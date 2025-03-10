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

from streamlit_extras.grid import grid
    

###############################
###### S T R E A M L I T ######
###############################

st. set_page_config(layout="wide")

# Set header title
st.title('Library Builder')

#Define the topologies that populate the node selection dropdown
DefaultNodes=['ABC', 'AbC', '!ABC!', 'Ab(C)', 'A!bCD!', 'A!b(C)DE!', 'Specify Other' ]

#These are the values that populate the dropdowns in col2 and end up ingested as part of the `BuildingBlock.FG` attribute.   
reactivity_types =   ['none',
                    'alkyl_amine',
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
                    'carboxylic_acid_ester',
                    'boronic_acid',
                    'boronic_ester',
                    'vinyl_alkene',
                    'azide',]

reaction_types = ['amide_coupling',
                  'amine_displacement',
                  'alcohol_displacement',
                  'Suzuki_coupling',
                  'reductive amination',
                  'heck_coupling',
                  'buchwald',
                  'wittig_olefination',
                  'sonogashira_coupling',
                  'click_reaction']
                                
reactivity_dict = {'none':['none'],
                  'amide_coupling': ['alkyl_amine', 'aryl_amine', 'carboxylic_acid'],
                  'fmoc_depro-amide_coupling': ['fmoc_amine', 'carboxylic_acid'],
                  'ester_depro-inv_amide_coupling': ['carboxylic_acid_ester', 'alkyl_amine', 'aryl_amine',],
                  'amine_displacement': ['alkyl_amine', 'aryl_amine','alkyl_halide'],
                  'alcohol_displacement': ['alkyl_alcohol', 'aryl_alcohol', 'alkyl_halide'],
                  'Suzuki_coupling' : ['boronic_acid', 'aryl_halide'],
                  'reductive amination': ['alkyl_amine', 'aryl_amine', 'fmoc_amine', 'aldehyde'],
                  'heck_coupling': ['vinyl_alkene', 'aryl_halide'],
                  'buchwald': ['alkyl_amine', 'aryl_amine', 'fmoc_amine', 'aryl_halide'],
                  'wittig_olefination': ['aldehyde', 'ketone', 'alkyl_halide'],
                  'sonogashira_coupling': ['aryl_halide', 'alkyl_halide', 'alkyne'],
                  'click_reaction':['azide', 'alkyne']}


class BuildingBlock:
    def __init__(self, BBid:str, valence:int, FGs:list, rxn:list):
        self.id = BBid
        self.label = BBid
        self.valence = valence
        self.FGs = FGs
        
        self.valence_nodes = [f"{chr(65+i)}{self.id[-1]}" for i in range(valence)]
        self.valence_dict = dict(zip(self.valence_nodes, self.FGs))
        self.valence_coords = []
        if len(rxn) != valence:
            self.rxn = rxn * valence
        else:
            self.rxn = rxn


        self.ValenceMap = {}
    
        self.calc_bb_coords(400)
    
    def add_bb_edge(self, node1:str, node2:str):
           self.bb_edges.append [node1, node2]
    def add_valence_nodes(self):
        for i in range(self.valence):
            bb_nodes.append (f"{chr(65+i)}{BB.id[-1]}")
    

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


class DEL_Graph:
    def __init__(self, n_BBs):
        self.Bead = BuildingBlock('Bead')
        self.BB_dict = {'Bead': self.Bead}

        try:
            n_BBs = int(n_BBs)  # Convert to integer if it's a string
        except ValueError:
            raise ValueError("n_BBs must be an integer or a string representing an integer")

        for i in range(1, n_BBs + 1):
            self.BB_dict[f'BB{i}'] = BuildingBlock(f'BB{i}')

def get_possible_rxns(reactant:str, reactivity_dict:dict):
    possible_reactions = []
    for reaction, reactants in reactivity_dict.items():
        if reactant in reactants:
            possible_reactions.append(reaction)
    return(possible_reactions)



#####################
## SIDEBAR OPTIONS ##
#####################

# Number of BB rows to display
num_rows = st.sidebar.number_input("Number of BBs:", min_value=1, max_value=10, value=3)

#topology = st.sidebar.radio("What topology?",["Linear", "Branched"],)

#bead_reactivity = st.sidebar.radio("What bead reactivity?",["alkyl_amine", "carboxylic_acid"],)




####################
## COLUMN OPTIONS ##
####################



# Store user inputs in a dict
data = {}

#Use separate handling for the Bead (we assume the bead has a valence of 1 and only uses amine or carboxyl FGs)
st.write(f"#### Bead Options ")  

# Create columns for different input types
col1, col2, col3, col4 = st.columns([1,1,2,1])

bead_valence = col1.number_input(f"", min_value=1, max_value=1, value=1, key=f"valence_bead")
bead_input = col2.selectbox(f"Bead handle ", ['alkyl_amine', 'carboxylic_acid'], key=f"fg_Bead") 
bead_reactions = get_possible_rxns(bead_input, reactivity_dict)
bead_reactions = col3.pills("Possible Reactions:", bead_reactions, selection_mode='multi', key = f'rxn_bead') 
data['Bead'] = {'valence':bead_valence, 'fg': bead_input, 'rxn': bead_reactions,}
for i in range(num_rows):
    
    
    
    FG_inputs = []

    
    st.write(f"#### BB{i+1} Options ")  # Section for each row
    # Create columns for different input types
    col1, col2, col3, col4 = st.columns([1,1,2,1])
    

    # Input for BB Valence
    bb_valence = col1.number_input(f"Valence (BB{i+1})", min_value=1, max_value=3, value=1, key=f"valence_{i}")  # Text input
    

    #Input for BB Functional Groups (FGs)
    for j in range(int(bb_valence)):
        FG_input = col2.selectbox(f"Reactive handle {j+1}", reactivity_types, key=f"fg_{i}.{j}")  # Generate input fields dynamically
        FG_inputs.append(FG_input)

    #Input for reactions
    possible_reactions = []
    for fg in FG_inputs:
        reactions = get_possible_rxns(fg, reactivity_dict)
        possible_reactions.extend(reactions)
    possible_reactions = list(set(possible_reactions))
    selected_reactions = col3.pills("Possible Reactions:", possible_reactions, selection_mode='multi', key = f'rxn_{i}') 
    #selected_reactions = col3.multiselect("Possible Reactions:",  possible_reactions, placeholder="Choose from the available reactions", key = f'rxn_{i}')
    
    # Append row data
    data[f'BB{i+1}'] = {'valence':bb_valence, 'fg': FG_inputs, 'rxn': selected_reactions,}


    with col4:
        pass


import streamlit as st
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
import streamlit.components.v1 as components



### -------------------- NETWORK GRAPH GENERATION -------------------- ###
col1, col2, col3, col4 = st.columns([2,1,1,2])
# Create a NetworkX graph
G = nx.Graph()

# Add nodes (BBs)
for bb, details in data.items():
    G.add_node(bb, label=bb)

# Add edges based on common reactions
for bb1, details1 in data.items():
    for bb2, details2 in data.items():
        if bb1 != bb2:  # Avoid self-loops
            common_reactions = set(details1["rxn"]) & set(details2["rxn"])
            if common_reactions:  # If they share at least one reaction
                G.add_edge(bb1, bb2, label=", ".join(common_reactions))


with col1:
    with st.container(border=True):
        # Optional: Display Matplotlib Graph (Static)
        st.subheader("BB Network Graph (Static)")
        plt.figure(figsize=(4, 3))
        pos = nx.spring_layout(G)  # Positioning of nodes
        nx.draw(G, pos, with_labels=True, node_color="skyblue", edge_color="gray", node_size=1000, font_size=10)

        # Draw edge labels (reaction names)
        edge_labels = {(u, v): d["label"] for u, v, d in G.edges(data=True)}
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=6)

        # Show the plot in Streamlit
        st.pyplot(plt)

with col2:
    # Display the entered data
    st.subheader("User Entries:")
    st.write(data)

with col3:
    # Generate adjacency matrix as a NumPy array
    adj_matrix = nx.to_numpy_array(G, dtype=int)

    # Convert to Pandas DataFrame for better visualization
    nodes = list(G.nodes)
    adjacency_df = pd.DataFrame(adj_matrix, index=nodes, columns=nodes)
    with st.expander("View Adjacency Matrix"):
        st.dataframe(adjacency_df)
with col4:
    pass
