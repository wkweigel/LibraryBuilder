
#####################################
###### D E P E N D A N C I E S ######
#####################################

import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import networkx as nx
import matplotlib.pylab as plt
import numpy as np
from pyvis.network import Network
from PIL import Image
from collections import defaultdict
from queue import Queue
import string
import math


    

###############################
###### S T R E A M L I T ######
###############################

st. set_page_config(layout="wide")

# Set header title
st.title('Library Builder')

#Define the topologies that populate the node selection dropdown
DefaultNodes=['ABC', 'AbC', '!ABC!', 'Ab(C)', 'A!bCD!', 'A!b(C)DE!', 'Specify Other' ]

reactivity_types =   ['alkyl_amine',
                    'aryl_amine',
                    'boc_amine',
                    'fmoc_amine',
                    'alkyl_alcohol',
                    'aryl_alcohol',
                    'alkyl_halide',
                    'aryl_halide',
                    'aldehyde',
                    'ketone',
                    'alkyne',
                    'carboxylic_acid',
                    'boronic_acid',
                    'boronic_ester',
                    'vinyl_alkene',
                    'azide',]

reaction_types = ['amide_bond_formation',
                  'amine_displacement',
                  'alcohol_displacement',
                  'Suzuki_coupling',
                  'reductive amination',
                  'heck_coupling',
                  'buchwald_hartwig_coupling',
                  'wittig_olefination',
                  'sonogashira_coupling',
                  'click_reaction']
                                
reactivity_dict = {'amide_bond_formation': ['alkyl_amine', 'aryl_amine','carboxylic_acid'],
                  'amine_displacement': ['alkyl_amine', 'aryl_amine','alkyl_halide'],
                  'alcohol_displacement': ['alkyl_alcohol', 'aryl_alcohol', 'alkyl_halide'],
                  'Suzuki_coupling' : ['boronic_acid', 'aryl_halide'],
                  'reductive amination': ['alkyl_amine', 'aryl_amine','aldehyde'],
                  'heck_coupling': ['vinyl_alkene', 'aryl_halide'],
                  'buchwald_hartwig_coupling': ['alkyl_amine', 'aryl_amine','aryl_halide'],
                  'wittig_olefination': ['aldehyde', 'ketone', 'alkyl_halide'],
                  'sonogashira_coupling': ['aryl_halide', 'alkyl_halide', 'alkyne'],
                  'click_reaction':['azide', 'alkyne']}

def get_possible_rxns(reactant:str, reactivity_dict:dict):
    possible_reactions = []
    for reaction, reactants in reactivity_dict.items():
        if reactant in reactants:
            possible_reactions.append(reaction)
    return(possible_reactions)

def get_possible_reactants(reaction:str, reactivity_dict:dict, chosen_reagent:str = None):
    possible_reagents = reactivity_dict[reaction]
    if chosen_reagent != None:
        possible_reagents = [reagent for reagent in possible_reagents if reagent != chosen_reagent]
    return possible_reagents

class BuildingBlock:
    def __init__(self, BBid:str):
        self.id=BBid
        self.label=BBid
        self.valence=None
        self.valence_coords = []
    
        self.reactivity_A=None
        self.reactivity_B=None
        self.reactivity_C=None

        self.ValenceMap = {}
    
        self.calc_bb_coords(400)
    

    def calc_bb_coords(self, multiplier:int):
        if self.id == 'Bead':
            self.bb_x = 0.0
            self.bb_y = 0.0
        else:
            self.bb_x = (int(self.id[-1]) * multiplier)
            self.bb_y = 0.0

    def calculate_valence_coords(self, cx, cy, n, r):
        """
        Calculate n equidistant points around a center (cx, cy) at radius r.
        
        Args:
        cx (float): X-coordinate of center
        cy (float): Y-coordinate of center
        n (int): Number of points
        r (float): Radius (distance from center)
        
        Returns:
        list of lists: [[x1, y1], [x2, y2], ..., [xn, yn]]
        """
        points = []
        for i in range(n):
            theta = (2 * math.pi * i) / n  # Evenly spaced angles
            x = cx + r * math.cos(theta)
            y = cy + r * math.sin(theta)
            self.valence_coords.append([x, y])
        
    


    def add_reactivities_list(self, reactivity_list):
        if self.BBvalence is None:
            return 'Valence Error'
        if self.valence == 'monovalent' and self.reactivity_A != None:
            self.reactivity_A = reactivity_list[0]
        if self.valence == 'divalent' and len(reactivity_list) == 2:
            self.reactivity_A = reactivity_list[0]
            self.reactivity_B = reactivity_list[1]
        if self.valence == 'trivalent' and len(reactivity_list) == 3:
            self.reactivity_A = reactivity_list[0]
            self.reactivity_B = reactivity_list[1]
            self.reactivity_C = reactivity_list[2]

    def add_reactivities(self, reactivity_dict):
        if self.valence is None:
            return 'Valence Error'
        if self.valence == 'monovalent' and len(reactivity_dict) == 1:
            self.reactivity_A = reactivity_dict['A']
        if self.valence == 'divalent' and len(reactivity_dict) == 2:
            self.reactivity_A = reactivity_dict['A']
            self.reactivity_B = reactivity_dict['B']
        if self.valence == 'trivalent' and len(reactivity_dict) == 3:
            self.reactivity_A = reactivity_dict['A']
            self.reactivity_B = reactivity_dict['B']
            self.reactivity_C = reactivity_dict['C']

    



class DEL_Graph:
    def __init__(self, n_BBs):
        self.Bead = BuildingBlock('Bead')
        if n_BBs == '1':
            self.BB1 = BuildingBlock('BB1')
            self.BB_dict = {'Bead': self.Bead, 'BB1' : self.BB1} 
        if n_BBs == '2':
            self.BB1 = BuildingBlock('BB1')
            self.BB2 = BuildingBlock('BB2')
            self.BB_dict = {'Bead': self.Bead, 'BB1' : self.BB1, 'BB2' : self.BB2} 
        if n_BBs == '3':
            self.BB1 = BuildingBlock('BB1')
            self.BB2 = BuildingBlock('BB2')
            self.BB3 = BuildingBlock('BB3')
            self.BB_dict = {'Bead': self.Bead, 'BB1' : self.BB1, 'BB2' : self.BB2, 'BB3' : self.BB3} 
        if n_BBs == '4':
            self.BB1 = BuildingBlock('BB1')
            self.BB2 = BuildingBlock('BB2')
            self.BB3 = BuildingBlock('BB3')
            self.BB4 = BuildingBlock('BB4')
            self.BB_dict = {'Bead': self.Bead, 'BB1' : self.BB1, 'BB2' : self.BB2, 'BB3' : self.BB3, 'BB4' : self.BB4} 
        
        self.Nodes = ['Bead'] + [f'BB{i}' for i in range(1, int(n_BBs) + 1)]
        self.building_blocks = {name: BuildingBlock(name) for name in self.Nodes}

        self.Edges = []
        self.Reactions = {}
        self.Labels = self.Nodes
        self.BB1_valence = None
        self.BB2_valence = None
        self.BB3_valence = None
        self.BB4_valence = None
        #self.Make_Graph()
        self.Update_Graph()

    
    
    def Make_Graph(self):
        '''
        Creates lists for the colors and shapes for each element node in the tree
        '''
        valence_dict = {'A': 0, 'B': 1, 'C': 2}
        self.node_colors=[]
        self.node_shapes=[]
        self.node_sizes = []
        #self.Labels = self.Nodes
        self.node_x = []
        self.node_y = []
        for node in self.Nodes:
            if node == 'Bead':
                self.node_colors.append('#3A3B3C')
                self.node_shapes.append('dot')
                self.node_sizes.append(400)
                self.node_x.append(0.0)
                self.node_y.append(0.0)

            if len(node) == 2:
                self.node_colors.append('#012A4A')
                self.node_shapes.append('hexagon')
                self.node_sizes.append(30)
                if node[-1] == '0':
                    self.node_x.append(50.0)
                    self.node_y.append(0.0)
                if node[-1] == '1': # This is for a BB1 Node
                    valence_coords = self.BB1.valence_coords #These are the valence coords associated with BB1
                    coord = valence_coords[valence_dict[node[0]]] #This the coord in the list of coords asscociated with the current valence node
                    self.node_x.append(coord[0])
                    self.node_y.append(coord[1])
                if node[-1] == '2': # This is for a BB2 Node
                    valence_coords = self.BB2.valence_coords #These are the valence coords associated with BB2
                    coord = valence_coords[valence_dict[node[0]]] #This the coord in the list of coords asscociated with the current valence node
                    self.node_x.append(coord[0])
                    self.node_y.append(coord[1])
                if node[-1] == '3': # This is for a BB3 Node
                    valence_coords = self.BB3.valence_coords #These are the valence coords associated with BB3
                    coord = valence_coords[valence_dict[node[0]]] #This the coord in the list of coords asscociated with the current valence node
                    self.node_x.append(coord[0])
                    self.node_y.append(coord[1])
                if node[-1] == '4': # This is for a BB4 Node
                    valence_coords = self.BB4.valence_coords #These are the valence coords associated with BB4
                    coord = valence_coords[valence_dict[node[0]]] #This the coord in the list of coords asscociated with the current valence node
                    self.node_x.append(coord[0])
                    self.node_y.append(coord[1])
            
            elif node in reaction_types:
                self.node_colors.append('#FA8072')
                self.node_shapes.append('square')
                self.node_sizes.append(200)


                
            elif node in list(self.BB_dict.keys()):
                if node[-1] == '1':
                    self.node_colors.append('#A9D6E5')
                    self.node_shapes.append('dot')
                    self.node_sizes.append(300)
                    self.node_x.append(self.BB1.bb_x)
                    self.node_y.append(self.BB1.bb_y)
                if node[-1] == '2':
                    self.node_colors.append('#61A5C2')
                    self.node_shapes.append('dot')
                    self.node_sizes.append(300)
                    self.node_x.append(self.BB2.bb_x)
                    self.node_y.append(self.BB2.bb_y)
                if node[-1] == '3':
                    self.node_colors.append('#468FAF')
                    self.node_shapes.append('dot')
                    self.node_sizes.append(300)
                    self.node_x.append(self.BB3.bb_x)
                    self.node_y.append(self.BB3.bb_y)
                if node[-1] == '4':
                    self.node_colors.append('#2A6F97')
                    self.node_shapes.append('dot')
                    self.node_sizes.append(300)
                    self.node_x.append(self.BB4.bb_x)
                    self.node_y.append(self.BB4.bb_y)

    def Replace_Node(self, current_node:str, replacement_node:str):
        self.Nodes = [replacement_node if node == current_node else node for node in self.Nodes]
        for i, edge in enumerate(self.Edges):
            for j, node in enumerate(edge):
                if node == current_node:
                    self.Edges[i][j] = replacement_node
        self.Update_Graph()

    def Replace_Label(self, current_node:str, replacement_label:str):
        self.Labels = [replacement_label if node == current_node else node for node in self.Labels]
        self.Update_Graph()

    def Update_Graph(self):
        #Initialize Graph object
        self.Graph= Network()

        #load the nodes, shapes, colors, and sizes for the DEL into the neworkx object
        self.Make_Graph()
        self.Graph.add_nodes(self.Nodes, 
                             color=self.node_colors, 
                             shape=self.node_shapes, 
                             value=self.node_sizes, 
                             label=self.Labels, 
                             #x=self.node_x, 
                             #y=self.node_y
                             )

        #load the edges for the DEL into the neworkx object
        self.Graph.add_edges(self.Edges)


def add_valence_nodes(Graph: DEL_Graph, BB: BuildingBlock, valence_type: str):
    '''Add the proper number of valence nodes to a BuildingBlock based on its specified valence type.'''
    BB_digit = BB.id[-1]
    
    # Extract the valence number from the string using a mapping dict
    valence_map = {"monovalent": 1, "divalent": 2, "trivalent": 3}
    valence_key = valence_type.replace(f"{BB.id}-", "")

    if valence_key in valence_map:

        BB.valence = valence_key
        BB.calculate_valence_coords(BB.bb_x, BB.bb_y, valence_map[valence_key], 50)

        for i in range(valence_map[valence_key]):
            node_id = f"{chr(65 + i)}{BB_digit}"  # Generates 'A1', 'B1', 'C1', etc. by enumerating starting with character "A"
            Graph.Nodes.append(node_id)
            Graph.Labels.append(node_id)
            Graph.Edges.append([f"BB{BB_digit}", node_id])
    elif valence_type == 'Bead':
        node_id = 'A0'
        Graph.Nodes.append(node_id)
        Graph.Edges.append(["Bead", node_id])

    return (Graph.Nodes, Graph.Edges)

def show_valence_options(Graph: DEL_Graph, BB: BuildingBlock, valence_type: str):
    """
    Displays Streamlit select boxes for choosing reactivity options based on the valence type 
    of a given BuildingBlock (BB). Updates the apropriate BulidingBlock class with intended reactivities. Creates the neccessary  
    valence nodes and updates the nodes and edges for the provided Graph.

    Parameters:
    -----------
    Graph : DEL_Graph
        The graph where nodes and edges representing valence interactions will be added.
    BB : BuildingBlock
        The building block whose valence and reactivities are being defined.
    valence_type : str
        The type of valence, formatted as "{BB.id}-monovalent", "{BB.id}-divalent", or "{BB.id}-trivalent".
    """
    valence_map = {"monovalent": 1, "divalent": 2, "trivalent": 3}

    # Extract the valence type (e.g., 'monovalent', 'divalent', 'trivalent')
    valence_key = valence_type.replace(f"{BB.id}-", "")
    
    if valence_key in valence_map:
        num_reactivities = valence_map[valence_key]

        # Create select boxes based on valence count
        reactivities = [st.selectbox(f'{BB.id} Reactivity {chr(65+i)}:', ['choose...'] + reactivity_types) for i in range(num_reactivities)]

        # Assign reactivities
        BB.reactivity_A, BB.reactivity_B, BB.reactivity_C = (reactivities + [None] * (3 - num_reactivities))[:3]

        # Update ValenceMap 
        for i, reactivity in enumerate(reactivities):
            BB.ValenceMap[reactivity] = f"{chr(65+i)}{BB.id[-1]}"

        # Add nodes and edges 
        Graph.Nodes, Graph.Edges = add_valence_nodes(Graph, BB, valence_type)

    return Graph.Nodes, Graph.Edges


# Function to add BB and reaction selections
def add_reaction_options():
    if len(st.session_state["pill_BBs"]) == 2:
        bb_A = DEL.BB_dict[st.session_state["pill_BBs"][0]]
        bb_B = DEL.BB_dict[st.session_state["pill_BBs"][1]]
        
        #Get lists for the two reactive handles on each of the BBs
        reactivityA_options = get_BB_reactivities(bb_A)
        reactivityB_options = get_BB_reactivities(bb_B)
        
        #Get lists of the reactions each of the reactive handles could participate in
        reactionsA = [rxn for reactivity in reactivityA_options for rxn in get_possible_rxns(reactivity, reactivity_dict)]
        reactionsB = [rxn for reactivity in reactivityB_options for rxn in get_possible_rxns(reactivity, reactivity_dict)]
        
        #Get the common reaction(s) between both reaction lists
        possible_reactions = list(set(reactionsA) & set(reactionsB))

        
        if possible_reactions:
            # If multiple reactions exist, use selectbox to choose one
            selected_rxn = st.selectbox("Select a reaction:", possible_reactions) if len(possible_reactions) > 1 else possible_reactions[0]
            

            # Append to session state lists
            st.session_state["bb_items"].append(st.session_state["pill_BBs"])
            st.session_state["rxn_items"].append(selected_rxn)

# Function to determine possible reactions for a reactant
def get_possible_rxns(reactant: str, reactivity_dict: dict):
    return [reaction for reaction, reactants in reactivity_dict.items() if reactant in reactants]

def get_BB_reactivities(BB):
    reactivities = []
    if BB.valence == 'monovalent':
        reactivities = [BB.reactivity_A]
    elif BB.valence == 'divalent':
        reactivities = [BB.reactivity_A, BB.reactivity_B]
    elif BB.valence == 'trivalent':
        reactivities = [BB.reactivity_A, BB.reactivity_B, BB.reactivity_C]
    return reactivities



#####################
## SIDEBAR OPTIONS ##
#####################

n_BBs = st.sidebar.radio(
"How many building blocks?",
["2", "3", "4"],)

topology = st.sidebar.radio(
"What topology?",
["Linear", "Branched"],)

bead_reactivity = st.sidebar.radio(
"What bead reactivity?",
["alkyl_amine", "carboxylic_acid"],)


#Define the nodes and edges for initial graph construction 
#These values are then updated based on user input 
Nodes = ['Bead', 'BB1', 'BB2']
Edges = []

BB3_valence= 'choose...'
BB4_valence= 'choose...'


DEL=DEL_Graph(n_BBs)
Edges = []
DEL.Bead.valence = 'monovalent'
DEL.Bead.reactivity_A = bead_reactivity
DEL.Bead.ValenceMap[bead_reactivity] = 'A0'
DEL.Nodes, DEL.Edges = add_valence_nodes(DEL, DEL.Bead, 'Bead')
DEL.Replace_Label('A0', DEL.Bead.reactivity_A)


####################
## COLUMN OPTIONS ##
####################

# Initialize session state for storing selections if not set
if "bb_items" not in st.session_state:
    st.session_state["bb_items"] = []
if "rxn_items" not in st.session_state:
    st.session_state["rxn_items"] = []
if "pill_BBs" not in st.session_state:
    st.session_state["pill_BBs"] = []  # Track selected BBs with pills input

# The columns for horizontal input
col1, col2, col3, col4 = st.columns([1, 1, 1, 3])

with col1:
    with st.container():
        st.subheader("Specify Valencies...")
        if n_BBs=='2':
            BB1_valence= st.selectbox('BB1 Valence:', ['choose...', 'BB1-monovalent', 'BB1-divalent', 'BB1-trivalent'])
            BB2_valence= st.selectbox('BB2 Valence:', ['choose...', 'BB2-monovalent', 'BB2-divalent', 'BB2-trivalent'])
        if n_BBs=='3':
            BB1_valence= st.selectbox('BB1 Valence:', ['choose...', 'BB1-monovalent', 'BB1-divalent', 'BB1-trivalent'])
            BB2_valence= st.selectbox('BB2 Valence:', ['choose...', 'BB2-monovalent', 'BB2-divalent', 'BB2-trivalent'])
            BB3_valence= st.selectbox('BB3 Valence:', ['choose...', 'BB3-monovalent', 'BB3-divalent', 'BB3-trivalent'])
        if n_BBs=='4':
            BB1_valence= st.selectbox('BB1 Valence:', ['choose...', 'BB1-monovalent', 'BB1-divalent', 'BB1-trivalent'])
            BB2_valence= st.selectbox('BB2 Valence:', ['choose...', 'BB2-monovalent', 'BB2-divalent', 'BB2-trivalent'])
            BB3_valence= st.selectbox('BB3 Valence:', ['choose...', 'BB3-monovalent', 'BB3-divalent', 'BB3-trivalent'])
            BB4_valence= st.selectbox('BB4 Valence:', ['choose...', 'BB4-monovalent', 'BB4-divalent', 'BB4-trivalent'])

with col2:
    st.subheader("Specify Reactive Handles...")
    if BB1_valence != 'choose...':
        st.subheader("BB1 Reactivities:")
        DEL.Nodes, DEL.Edges = show_valence_options(DEL, DEL.BB1, BB1_valence)
    if BB2_valence != 'choose...':
        st.subheader("BB2 Reactivities:")
        DEL.Nodes, DEL.Edges = show_valence_options(DEL, DEL.BB2, BB2_valence)
    if BB3_valence != 'choose...':
        st.subheader("BB3 Reactivities:")
        DEL.Nodes, DEL.Edges = show_valence_options(DEL, DEL.BB3, BB3_valence)
    if BB4_valence != 'choose...':
        st.subheader("BB4 Reactivities:")
        DEL.Nodes, DEL.Edges = show_valence_options(DEL, DEL.BB4, BB4_valence)




with col3:
    st.subheader("Specify Reaction Participants...")
    # Use st.pills for selection, ensuring max 2 selections
    BB_options = list(DEL.BB_dict.keys())
    st.session_state["pill_BBs"] = st.pills("Reaction Participants (Choose only two at a time):", BB_options, selection_mode='multi')  
    #Use a button to submit BB choices 
    if len(st.session_state["pill_BBs"]) == 2:
        st.button("Add BBs", on_click=add_reaction_options, type='primary')
    else:
        if st.button("Add BBs"):
            st.write('Input requires two BBs.')
        




# Layout
#with st.container():
    #st.button("Add BBs", on_click=add_reaction_options)

col_a, col_b, col_c = st.columns([0.3, 0.3, 3])

    
with col_a:
    st.subheader("BB Selections")
    for item in st.session_state["bb_items"]:
        st.write(", ".join(item))  # Display added BB selections

with col_b:
    st.subheader("Rxn Selections")
    for item in st.session_state["rxn_items"]:
        st.write(item)  # Display selected reactions

with col_c:
    pass


#####################
## DYNAMIC UPDATES ##
#####################

for bb_pair, rxn in zip(st.session_state["bb_items"], st.session_state["rxn_items"]):

    bb_A = DEL.BB_dict[bb_pair[0]]
    bb_B = DEL.BB_dict[bb_pair[1]]
    
    #Get lists for the two reactive handles on each of the BBs
    reactivityA_options = get_BB_reactivities(bb_A)
    reactivityB_options = get_BB_reactivities(bb_B)
    
    #Get lists of the reactions each of the reactive handles could participate in
    reactionsA = [rxn for reactivity in reactivityA_options for rxn in get_possible_rxns(reactivity, reactivity_dict)]
    reactionsB = [rxn for reactivity in reactivityB_options for rxn in get_possible_rxns(reactivity, reactivity_dict)]

    DEL.Reactions[rxn] = st.session_state["pill_BBs"]
    DEL.Nodes.append(rxn)
    DEL.Labels.append(rxn)
    target_reactivities = reactivity_dict[rxn]
    possible_reactions = list(set(reactionsA) & set(reactionsB))

    reactivityA_node=list(set(reactivityA_options) & set(target_reactivities))
    reactivityB_node=list(set(reactivityB_options) & set(target_reactivities))

    if len(reactivityA_node) == 1:
        DEL.Edges.append([rxn, bb_A.ValenceMap[reactivityA_node[0]]])
    if len(reactivityB_node) == 1:
        DEL.Edges.append([rxn, bb_B.ValenceMap[reactivityB_node[0]]])
        

# Replace the default valence node labels once their reactivity type has been chosen from the dropdown input
if int(n_BBs) >= 2:
    if DEL.BB1.reactivity_A != None:
        DEL.Replace_Label('A1', DEL.BB1.reactivity_A)
    if DEL.BB1.reactivity_B != None:
        DEL.Replace_Label('B1', DEL.BB1.reactivity_B)
    if DEL.BB1.reactivity_C != None:
        DEL.Replace_Label('C1', DEL.BB1.reactivity_C)

    if DEL.BB2.reactivity_A != None:
        DEL.Replace_Label('A2', DEL.BB2.reactivity_A)
    if DEL.BB2.reactivity_B != None:
        DEL.Replace_Label('B2', DEL.BB2.reactivity_B)
    if DEL.BB2.reactivity_C != None:
        DEL.Replace_Label('C2', DEL.BB2.reactivity_C)

if int(n_BBs) >= 3:
    if DEL.BB3.reactivity_A != None:
        DEL.Replace_Label('A3', DEL.BB3.reactivity_A)
    if DEL.BB3.reactivity_B != None:
        DEL.Replace_Label('B3', DEL.BB3.reactivity_B)
    if DEL.BB3.reactivity_C != None:
        DEL.Replace_Label('C3', DEL.BB3.reactivity_C)

if int(n_BBs) >= 4:
    if DEL.BB4.reactivity_A != None:
        DEL.Replace_Label('A4', DEL.BB4.reactivity_A)
    if DEL.BB4.reactivity_B != None:
        DEL.Replace_Label('B4', DEL.BB4.reactivity_B)
    if DEL.BB4.reactivity_C != None:
        DEL.Replace_Label('C4', DEL.BB4.reactivity_C)


DEL.Update_Graph()


#DEL.Graph.show_buttons(filter_=['physics'])


def display_pyvis_graph(network: Network):
    path = './'
    network.save_graph(f'{path}/del.html')
    st.components.v1.html(open(f'{path}/del.html', "r").read(), height=600, width=1000, scrolling=False)

with col4:
    st.subheader("Library Viewer")
    display_pyvis_graph(DEL.Graph)



def generate_adjacency_matrix(pyvis_net: Network):
    """
    Converts a Pyvis network graph into an adjacency matrix.

    Parameters:
    -----------
    pyvis_net : pyvis.network.Network
        The Pyvis network graph to be converted.

    Returns:
    --------
    adjacency_df : pandas.DataFrame
        A DataFrame representing the adjacency matrix.
    """
    # Convert Pyvis Network to a NetworkX graph (Pyvis uses NetworkX internally)
    nx_graph = nx.Graph()

    # Add nodes and edges from Pyvis graph
    for node in pyvis_net.nodes:
        nx_graph.add_node(node["id"])
    
    for edge in pyvis_net.edges:
        nx_graph.add_edge(edge["from"], edge["to"])

    # Generate adjacency matrix as a NumPy array
    adj_matrix = nx.to_numpy_array(nx_graph, dtype=int)

    # Convert to Pandas DataFrame for better visualization
    nodes = list(nx_graph.nodes)
    adjacency_df = pd.DataFrame(adj_matrix, index=nodes, columns=nodes)

    return adjacency_df

matrix = generate_adjacency_matrix(DEL.Graph)

matrix
