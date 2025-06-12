__version__ = "0.1.0"

# docopt CLI
__doc__ = """
Collaspi: Post-(Extracted) Layout Simulations Speed-Up through Netlist Reduction.

Usage:
    collaspi <netlist_file> <collapsed_netlist_file> [--maxr=MAX_RESISTANCE] [--maxc=MAX_CAPACITANCE]
    collaspi [--test]
    collaspi --version
    collaspi --help
Options:
    --test          Run the test suite.
    --version       Show the version of Collaspi.
    --help          Show this help message.
"""

from docopt import docopt
from tqdm import tqdm, trange

from pathlib import Path
from typing import Dict, Optional, Union, Tuple
from collections import defaultdict, UserDict
#from itertools import product #combinations
from networkx import (
    Graph, 
    resistance_distance,
    connected_components, 
    betweenness_centrality, 
    periphery
)

import PySpice.Logging.Logging as Logging
logger = Logging.setup_logging()
from PySpice.Spice.Netlist import Circuit, SubCircuit
from PySpice.Spice.Parser import SpiceParser
from PySpice.Spice.BasicElement import Resistor, Capacitor, BehavioralResistor, BehavioralCapacitor


class Suffix(UserDict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

_suffix = Suffix({
        't': 1e12, 'g': 1e9, 'meg': 1e6, 'k': 1e3,
        'm': 1e-3, 'u': 1e-6, 'n': 1e-9,
        'p': 1e-12, 'f': 1e-15, 'M': 1e6, 'K': 1e3,
        'G': 1e9, 'T': 1e12, 'P': 1e15
    })

class ReportConfig(object):
    maxr: Optional[float]
    maxc: Optional[float]
    unitr = 1.0
    unitc = _suffix['f']

def spice_to_float(val):
    if not isinstance(val, str):
        return float(val)
    
    suffixes = _suffix

    val = val.strip().lower()
    for suffix in sorted(suffixes, key=len, reverse=True):
        if val.endswith(suffix):
            return float(val[:-len(suffix)]) * suffixes[suffix]
    return float(val)

# Deprecated: this function is not required when building subcircuit from 
# Calibre PEX 'DSPF' extracted net file format.
def pre_process_netlist(path: Optional[Union[Path, str]]=None, source: Optional[str]=None) -> Union[Circuit, SubCircuit]:
    """Pre-process the netlist to ensure it is in a suitable format for analysis."""
    raise DeprecationWarning('This function is no longer required.')
    netlist: str = None
    if path is not None:
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Netlist file {path} does not exist.")
        with open(path, 'r') as f:
            netlist = f.read()
    elif source is not None:
        netlist = source
    else:
        raise ValueError("Either path or source must be provided.")
    
    #rename nodes to replace ':' with '_' - not necessary
    #netlist = netlist.replace(':', '_')
    
    # rename repeated element names to their node1_node2 format
    lines = netlist.splitlines()
    element_names = defaultdict(int)
    comments = []
    nets = []
    for i, line in tqdm(enumerate(lines), desc="Discovering Nets"):
        if line.startswith('**'):
            # detect net name
            net_name = line.split()[1]
            nets.append((net_name,i))
    current_net = nets.pop(0) if nets else None
    for i, line in tqdm(enumerate(lines), desc="Pre-Processing Netlist"):
        if line.startswith('.') or not line.strip():
            continue
        if line.startswith('+'):
            continue  # continuation line, skip it
        if line.startswith('**') and nets:  # assuming net names start with '**'
            current_net = nets.pop(0)
        if line.startswith('*'):  # assuming comments start with '*'
            comments.append(line[1:].strip())
            continue
        parts = line.split()
        if len(parts) < 3: continue
        node1, node2 = parts[1], parts[2]
        # check if node has numeric prefix
        prefix = node1.split(':')[0]
        
        if prefix.isnumeric() and prefix != '0':

            parts[1] = f"{current_net[0]}:{node1}"
        prefix = node2.split(':')[0]
        if prefix.isnumeric() and prefix != '0':
            parts[2] = f"{current_net[0]}:{node2}"
        lines[i] = ' '.join(parts)
        element_name = parts[0]
        element_names[element_name] += 1
        if element_names[element_name] > 1:
            new_name = f"{element_name}_{element_names[element_name]}"
            lines[i] = f"{new_name} {' '.join(parts[1:])}"
    netlist = '\n'.join(lines)
    # TESTING purposes only
    with open('./pre_processed_netlist.sp', 'w') as f:
        f.write(netlist)
    # parse the netlist
    cir = SpiceParser(source=netlist).build_circuit()
    #for comment in comments:
    #    cir.raw_spice += f"* {comment}\n"
    # collect comments from the netlist
    return cir

def is_resistor(elem):
    return isinstance(elem, Resistor) or isinstance(elem, BehavioralResistor)

def is_capacitor(elem):
    return isinstance(elem, Capacitor) or isinstance(elem, BehavioralCapacitor)

def build_rcc_graph(netlist: Union[Circuit, SubCircuit]) -> Tuple[Graph, Graph]:
    """Build a resistance-capcitance-coupling capacitance graph from a given circuit."""
    RG = Graph()
    CCG = Graph()
    for node in netlist.nodes:
        RG.add_node(node)
        CCG.add_node(node)
    
    for elem in tqdm(netlist.elements, desc="Building RCC Graphs from Netlist"):
        if is_resistor(elem) and not elem.has_parameter('model'):
            resistance = elem.resistance if elem.has_parameter('resistance') else elem.resistance_expression
                
            if RG.has_edge(elem.nodes[0], elem.nodes[1]):
                # compute parallel resistance if the edge already exists
                r = RG.edges[elem.nodes[0], elem.nodes[1]]['weight']
                val = spice_to_float(resistance)
                RG.edges[elem.nodes[0], elem.nodes[1]]['weight'] = r * val / (r + val)
            else:
                RG.add_edge(elem.nodes[0], elem.nodes[1], name=elem.name, weight=spice_to_float(resistance))
        
        if is_capacitor(elem) and not elem.has_parameter('model'):
            capacitance = elem.capacitance if elem.has_parameter('capacitance') else elem.capacitance_expression
            if CCG.has_edge(elem.nodes[0], elem.nodes[1]):
                # compute parallel capacitance if the edge already exists
                c = CCG.edges[elem.nodes[0], elem.nodes[1]]['weight']
                CCG.edges[elem.nodes[0], elem.nodes[1]]['weight'] = c + spice_to_float(capacitance)
            else:
                CCG.add_edge(elem.nodes[0], elem.nodes[1], name=elem.name, weight=spice_to_float(capacitance))
    return RG, CCG


def build_lumped_elements_graphs(RG: Graph, CCG: Graph) -> Dict[str, Tuple[Graph, Graph]]:

    nets = [RG.subgraph(c).copy() for c in connected_components(RG)]

    collapsed_nets = []
    for net in nets:
        collapsed_net = Graph()
        edge_nodes = periphery(net)
        betweeness_centrality_score = betweenness_centrality(net, normalized=True, endpoints=False)
        center_node = max(betweeness_centrality_score, key=betweeness_centrality_score.get)
        r_index = 0
        for edge_node in edge_nodes:
          if edge_node == center_node: continue
          collapsed_net.add_edge(edge_node, center_node, name=f'R{r_index}', weight=resistance_distance(net, center_node, edge_node))
          r_index += 1
        collapsed_nets.append((collapsed_net, center_node))

    # build capacitance node to net interface map
    cap_coupling_map = defaultdict(list)
    for edge in CCG.edges:
        node1 = edge[0]
        node2 = edge[1]
        net1 = [k for k, net in enumerate(nets) if node1 in net.nodes][0]
        net2 = [k for k, net in enumerate(nets) if node2 in net.nodes][0]
        net1_center_node = collapsed_nets[net1][1]
        net2_center_node = collapsed_nets[net2][1]
        cap_coupling_map[(net1_center_node, net2_center_node)].append(CCG.edges[edge]['weight'])
    
    reduced_resistance_nets = [reduced_net for reduced_net,_ in collapsed_nets]
    cap_coupling_map = {k: sum(v) for k, v in cap_coupling_map.items()}
    return reduced_resistance_nets, cap_coupling_map

def build_collapsed_netlist(
    netlist: Union[Circuit, SubCircuit], 
    CCG_lumped_map: Dict[Tuple, float], 
    RG_lumped_nets: Dict[str, Tuple[Graph, Graph]],
    cfg: Optional[ReportConfig] = None,
    rounding: Optional[int] = 4
    ) -> Union[SubCircuit, Circuit]:
    # create a copy of the analysed netlist
    collapsed_netlist = None
    new_name = None
    if hasattr(netlist, 'name'):
        new_name = f"{netlist.name}_collapsed"
    else:
        new_name = "None_collapsed"
    if isinstance(netlist, Circuit):
        collapsed_netlist = Circuit(new_name)
    elif isinstance(netlist, SubCircuit):
        collapsed_netlist = SubCircuit(new_name)
    else:
        raise TypeError("netlist must be a Circuit or SubCircuit instance")
    # count the total number of nodes in the circuit
    visited_cap = []
    r_index = 0
    
    critical_report = {}
    
    for net in tqdm(RG_lumped_nets, desc="Building Lumped Resistance Netlist"):
        for edge in net.edges:
            # introduce intermediate node
            r = net.edges[edge]['weight']
            collapsed_netlist.R(f'{r_index}', \
                edge[0].name, edge[1].name, round(r, rounding))
            if cfg:
                if hasattr(cfg, 'maxr'):
                    if r > cfg.maxr:
                        critical_report.update({
                            f'R{r_index}': (edge[0].name, edge[1].name, round(r, rounding))
                        })
            r_index += 1
    c_index = 0
    for cap in tqdm(CCG_lumped_map, desc="Building Lumped Capacitance Netlist"):
        if cap not in visited_cap:
            visited_cap.append(cap)
            collapsed_netlist.C(f'{c_index}', \
                cap[0].name, cap[1].name, CCG_lumped_map[cap])
            if cfg:
                if hasattr(cfg, 'maxc'):
                    if CCG_lumped_map[cap] > cfg.maxc:
                        critical_report.update({
                            f'C{c_index}': (cap[0].name, cap[1].name, CCG_lumped_map[cap])
                        })
            c_index += 1
    
    # copy the remaining elements in the original spice into the 
    # collapsed passive parasitics generated spice
    for elem in netlist.elements:
        if not is_resistor(elem) and not is_capacitor(elem):
            collapsed_netlist._add_element(elem)
    # copy the comments from the original netlist
    # FIXME: need to check PySpice documentation for this
    for line in str(netlist).splitlines(): 
        if line.startswith('*'):  # assuming comments start with '*'
            collapsed_netlist.raw_spice += line + '\n'
    # copy the subcircuits from the original netlist
    for subckt in netlist.subcircuits:
        collapsed_netlist.subcircuit(subckt)
    # copy the global nodes from the original netlist
    # copy the parameters from the original netlist
    return collapsed_netlist, critical_report
    
    
def test_build_rcc_graph():
    # Example test case for build_rcc_graph
    circuit = Circuit('Test Circuit')
    circuit.R('R1', 'n1', 'n2', 100)
    circuit.C('C1', 'n2', 'n3', 1e-6)
    RG, CCG = build_rcc_graph(circuit)
    assert len(RG.edges) == 1
    assert len(CCG.edges) == 1
        

def test_build_lumped_element_graphs():
    from pprint import pprint
    # Example test case for build_collapsed_resistance_graphs
    circuit = Circuit('Test Circuit')
    circuit.R('R1', 'n1', 'n2', '100u')
    circuit.R('R2', 'n2', 'n3', 200)
    circuit.R('R3', 'n1', 'n3', 300)  # Coupled resistance
    circuit.C('C1', 'n2', 'n3', 1e-12)
    circuit.C('C2', 'n1', 'n3', 2e-12)
    circuit.C('C3', 'n1', 0, 1e-12)  # Cap to ground
    circuit.C('C4', 'n2', 0, 1e-12)  # Cap to ground
    RG, CCG = build_rcc_graph(circuit)
    RG_sub_nets_map, CCG_lumped_map = build_lumped_elements_graphs(RG, CCG)
    assert len(RG_sub_nets_map[1].edges) == 2
    assert len(CCG_lumped_map) == 2  # C1, C2, C3, C4 should be lumped into 2 elements
    
    circuit = Circuit('Test Circuit 2')
    circuit.R('R1', 'n1', 'n2', 100)
    circuit.R('R2', 'n2', 'n3', 200)
    circuit.R('R3', 'n1', 'n3', 300)
    circuit.R('R4', 'n3', 'n4', 200)
    circuit.R('R5', 'n4', 'n5', 100)
    circuit.R('R6', 'n5', 'n6', 150) 
    # cap to ground
    circuit.C('C1', 'n2', 0, 1e-12)
    circuit.C('C2', 'n1', 0, 2e-12)
    circuit.C('C3', 'n4', 0, 1e-12)
    circuit.C('C4', 'n3', 0, 1e-12)
    
    RG, CCG = build_rcc_graph(circuit)
    RG_sub_nets_map, CCG_lumped_map = build_lumped_elements_graphs(RG, CCG)
    print("Resistance Graph:")
    pprint(RG_sub_nets_map[1].edges(data=True))
    print("Capacitance Coupling Graph:")
    pprint(CCG_lumped_map)
   
def test_build_collapsed_netlist():
    from PySpice.Spice.Netlist import SubCircuitFactory
    # Example test case for build_collapsed_netlist
    class PexNet(SubCircuitFactory):
        __name__ = 'test_extracted_net'
        NODES = ('n1', 'n9', 'n5', 'n6')
        def __init__(self):
            super().__init__()
            self.R(1, 'n1', 'n2', 1)
            self.R(2, 'n2', 'n3', 1)
            self.R(3, 'n3', 'n4', 1)
            self.R(4, 'n3', 'n5', 1)
            self.R(5, 'n6', 'n7', 1)
            self.R(6, 'n7', 'n8', 1)
            self.R(7, 'n8', 'n9', 1)

            # add coupled cap
            self.C(1, 'n1', 'n6',1)
            self.C(2, 'n2', 'n7',1)
            self.C(10, 'n1', 'n7',1)

            # add cap to ground
            # add cap to ground for net1
            self.C(3, 'n3', 0,1)
            self.C(4, 'n4', 0,1)
            self.C(5, 'n5', 0,1)
            # add cap to ground for net2
            self.C(6, 'n8', 0,1)
            self.C(7, 'n9', 0,1)
            self.C(8, 'n6', 0,1)
            self.C(9, 'n7', 0,1)
    circuit = PexNet()
    RG, CCG = build_rcc_graph(circuit)
    RG_sub_nets_map, CCG_lumped_map = build_lumped_elements_graphs(RG, CCG)
    collapsed_netlist = build_collapsed_netlist(circuit, CCG_lumped_map, RG_sub_nets_map)
    print(str(collapsed_netlist))
    assert isinstance(collapsed_netlist, SubCircuit)
    assert len(collapsed_netlist.elements) > 0
    
    correct_collapsed_netlist = """
.subckt None_collapsed
R0 n1 n3 2.0
R1 n3 n4 1.0
R2 n3 n5 1.0
R3 n9 n7 2.0
R4 n7 n6 1.0
C0 0 n3 3
C1 0 n7 4
C2 n3 n7 3
.ends None_collapsed
"""
    # FIXME: Need to find a way to compare the netlists properly by ignoring node names
    diff = set(str(collapsed_netlist).splitlines()) - set(correct_collapsed_netlist.splitlines())
    assert len(diff) == 0, f"Collapsed netlist does not match expected output. Differences: {diff}"

def test():
    test_build_rcc_graph()
    test_build_lumped_element_graphs()
    test_build_collapsed_netlist()
    print("All tests passed!")
    
def main():
    args = docopt(__doc__, version=f"Collaspi v{__version__}")
    if '--test' in args and args['--test']:
        print("Running tests...")
        test()
    if args['--version']:
        print(f"Collaspi version {__version__}")
        return
    input_file = Path(args['<netlist_file>']).resolve()
    output_file = Path(args['<collapsed_netlist_file>'])
    if not input_file.exists():
        raise FileNotFoundError(f"Input file {input_file} does not exist.")
    
    cfg = ReportConfig()
    if args.get('--maxc'):
        cfg.maxc = spice_to_float(args['--maxc'])
    if args.get('--maxr'):
        cfg.maxr = spice_to_float(args['--maxr'])
        
    #netlist = pre_process_netlist(path=input_file, source=None)
    netlist = list(SpiceParser(path=input_file).build_circuit().subcircuits)[0]
    RG, CCG = build_rcc_graph(netlist)
    # assert if RG is not empty
    # when RG is empty, the extracted netlist used C or C+CC options,
    # which means that all extracted capacitances are already lumped
    if RG.number_of_nodes() == 0:
        raise ValueError("The resistance graph is empty. C+CC PEX extraction already provides lumped capacitance (maximally reduced) netlists. Please check the input netlist to result from R, RC or RCC PEX extraction.")
    
    RG_sub_nets_map, CCG_lumped_map = build_lumped_elements_graphs(RG, CCG)
    collapsed_netlist, report = build_collapsed_netlist(netlist, CCG_lumped_map, RG_sub_nets_map, cfg=cfg)
    
    if report:
        output_dir = output_file.parent
        with open(output_dir / 'report.log', 'w') as f:
            f.write('Critical parasitics: \n')
            for node in report:
                f.write(f'{node}: {report[node]}\n')

    
    with open(output_file, 'w') as f:
        f.write(str(collapsed_netlist)+'.END')
    print(f"Collapsed netlist written to {output_file}")
    print()
    print("Note: For Cadence Vivado integration, remove '.title' statement in the generated lumped netlist.")
    
    print("Done! :)")
    # TODO: add support for reducing parallel and series MOSFET devices

if __name__ == "__main__":
    from pprint import pprint
    
    input_file = Path('./data/tv_dynamic_ls.pex.netlist').resolve()
    output_file = Path('./data/output.pex.netlist')
    if not input_file.exists():
        raise FileNotFoundError(f"Input file {input_file} does not exist.")
    
    cfg = ReportConfig()
        
    #netlist = subckt_to_ckt(path=input_file, source=None)
    netlist = list(SpiceParser(path=input_file).subcircuits)[0]
    netlist = netlist.build()
    RG, CCG = build_rcc_graph(netlist)
    # assert if RG is not empty
    # when RG is empty, the extracted netlist used C or C+CC options,
    # which means that all extracted capacitances are already lumped
    if RG.number_of_nodes() == 0:
        raise ValueError("The resistance graph is empty. C+CC PEX extraction already provides lumped capacitance (maximally reduced) netlists. Please check the input netlist to result from R, RC or RCC PEX extraction.")
    
    RG_sub_nets_map, CCG_lumped_map = build_lumped_elements_graphs(RG, CCG)
    collapsed_netlist, report = build_collapsed_netlist(netlist, CCG_lumped_map, RG_sub_nets_map, cfg=cfg)
    
    if report:
        output_dir = output_file.parent
        with open(output_dir / 'report.log', 'w') as f:
            f.write('Critical parasitics: \n')
            for node in report:
                f.write(f'{node}: {report[node]}\n')

    
    with open(output_file, 'w') as f:
        f.write('Post-Processed by Collaspi\n'+str(collapsed_netlist))
    print(f"Collapsed netlist written to {output_file}")
    print()
    print("Note: For Cadence Vivado integration, remove '.title' statement in the generated lumped netlist.")
    
    print("Done! :)")
    
    #pprint(str(pre_process_netlist('./data/testrcc.pex.netlist')))
    #test()

__all__ = [
    build_rcc_graph,
    build_lumped_elements_graphs,
    build_collapsed_netlist,
]