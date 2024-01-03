# IMPORT RELEVANT PACKAGES
import sys
from PyQt6.QtCore import *
from PyQt6.QtWidgets import *
from PyQt6.QtGui import *
import os
import shutil
import pandas as pd
import numpy as np
import seaborn as sns
import pickle
import sys
import multiprocessing as mp
import time

sys.path.append('./safepy/')

import safe

import matplotlib.pyplot as plt

import networkx as nx
import argparse


class PandasModel(QAbstractTableModel):
    def __init__(self, dataframe: pd.DataFrame, parent=None):
        QAbstractTableModel.__init__(self, parent)
        self._dataframe = dataframe

    def rowCount(self, parent=QModelIndex()) -> int:
        if parent == QModelIndex():
            return len(self._dataframe)

        return 0

    def columnCount(self, parent=QModelIndex()) -> int:
        if parent == QModelIndex():
            return len(self._dataframe.columns)
        return 0

    def data(self, index: QModelIndex, role=Qt.ItemDataRole):
        if not index.isValid():
            return None

        if role == Qt.ItemDataRole.DisplayRole:
            return str(self._dataframe.iloc[index.row(), index.column()])

        return None

    def headerData(
        self, section: int, orientation: Qt.Orientation, role: Qt.ItemDataRole
    ):
        if role == Qt.ItemDataRole.DisplayRole:
            if orientation == Qt.Orientation.Horizontal:
                return str(self._dataframe.columns[section])

            if orientation == Qt.Orientation.Vertical:
                return str(self._dataframe.index[section])

        return None
        
        

class Window(QMainWindow):
    def __init__(self):
        
        super().__init__(parent=None)
        
        wid = QWidget(self)
        self.setCentralWidget(wid)
        main_layout = QGridLayout()
        wid.setLayout(main_layout)
        
        self.setWindowTitle("ProteoML Genius Visualizer")   # label window
        
        
        
        # Setup constants for SAFE run
        self.go_matrix = []
        
        self.node_distance_metric = 'shortpath_weighted_layout'
        self.neighborhood_radius = 0.15
        
        self.multiple_testing = False
        self.background = 'network'
        self.full_force = True
        
        self.attribute_distance_threshold = 0.5
        
        self.background_color = '#000000'
        self.save_fig = "./runfiles/temp/"
        
        
        
        # Dropdown for network selection
        main_layout.addWidget(QLabel("Please select a base network for mapping:"), 0, 0, 1, 1)
        self.netselect = QComboBox()
        self.netselect.addItem("HuRI Human Proteome")   #0
        self.netselect.addItem("STRING-db Human Proteome")   #1
        main_layout.addWidget(self.netselect, 0, 1, 1, 1)
        
        #self.netselect.currentIndexChanged.connect(self.updateNetwork)   #user wants to change network view
        
        
        
        # Button for annotation upload
        self.autofit = QPushButton("Click to upload a file containing Uniprot IDs")
        main_layout.addWidget(self.autofit, 1, 0, 1, 2)
        
        main_layout.addWidget(QLabel(" "), 2, 0)   # spacer
        
        self.autofit.clicked.connect(self.upload)
        
        self.user_uids = None
        
        self.upload_notation = QLabel(" ")
        main_layout.addWidget(self.upload_notation, 3, 0)
        
        
        
        # Big box for displaying figure
        main_layout.addWidget(QLabel(" - Network Visualisation - "), 4, 0, 1, 1)
        self.sf_disp = QPixmap("./runfiles/vis_placeholder.tif")
        self.sf_label = QLabel()
        self.sf_disp_sm = self.sf_disp.scaledToWidth(800)
        self.sf_label.setPixmap(self.sf_disp_sm)
        main_layout.addWidget(self.sf_label, 5, 0, 1, 2)
        
        # Scrollable table for going through GO annotations
        self.go_lab = QLabel(" - GO Terms Mapped to Network - ")
        main_layout.addWidget(self.go_lab, 6, 0, 1, 1)
        
        self.mapped_go = pd.DataFrame([[0,0,'GO:0000000','None']], columns=['Proteins in Dataset', 'Proteins in Network', 'GO Term',  'Detailed Description'])
        self.go_model = PandasModel(self.mapped_go)
        self.view = QTableView()
        self.view.setModel(self.go_model)
        self.view.horizontalHeader().resizeSection(0, 150)
        self.view.horizontalHeader().resizeSection(1, 150)
        self.view.horizontalHeader().resizeSection(2, 150)
        self.view.horizontalHeader().resizeSection(3, 500)
        self.view.setMinimumSize(800, 150)
        main_layout.addWidget(self.view, 7, 0, 1, 2)
        
        
        
    def update_network(self):
        # User requested different network view - upload network
        if self.netselect.currentIndex() == 0:
            # HuRI view
            self.G = pickle.load(open('huri_network.pickle', 'rb'))
        elif self.netselect.currentIndex() == 1:
            # STRING-db view
            self.G = pickle.load(open('safe_network.pickle', 'rb'))
        # Display network as .pdf
        self.sf_label.setPixmap(QPixmap('path_to_network'))
        
    def upload(self):
        # User is providing an upload file with sites to map
        upload_from = QFileDialog.getOpenFileName(self, "Upload Data File Containing Uniprot IDs")
        upload_from_dir = upload_from[0]
            
        if ".xls" in upload_from_dir:
            # Excel file
            df = pd.read_excel(upload_from_dir)
            self.upload_notation.setText('Uploaded excel file: ' + upload_from_dir + '\nNow initiating SAFE run...')
        elif ".csv" in upload_from_dir:
            # csv file
            df = pd.read_csv(upload_from_dir)
            print('UPLOADED FILE')
            self.upload_notation.setText('Uploaded csv file: ' + upload_from_dir + '\nNow initiating SAFE run...')
        elif ".tsv" in upload_from_dir:
            df = pd.read_csv(upload_from_dir, sep='\t')
            self.upload_notation.setText('Uploaded tsv file: ' + upload_from_dir + '\nNow initiating SAFE run...')
        elif ".txt" in upload_from_dir:
            df = pd.read_csv(upload_from_dir, sep=" ")
            df = df.rename(columns={0:'ACC_ID'})
            self.upload_notation.setText('Uploaded txt file: ' + upload_from_dir + '\nNow initiating SAFE run...')
        else:
            # improper format >:(
            self.upload_notation.setText('Unable to upload file. Please attempt with a .xlsx, .csv, .txt, or .tsv file.\nCannot run SAFE until proper file format uploaded.')
            self.init_feat.setEnabled(False)
            df = None
        
        
        self.user_uids = df
        
        # Clean up data - pull just uniprot IDs
        if 'ACC_ID' in self.user_uids.columns:
            self.uids = self.user_uids['ACC_ID'].tolist()
        elif 'Uniprot ID' in self.user_uids.columns:
            self.uids = self.user_uids['Uniprot ID'].tolist()
            print(self.uids)
        elif 'uniprot_id' in self.user_uids.columns:
            self.uids = self.user_uids['uniprot_id'].tolist()
        else:
            # improper format >:(
            self.upload_notation.setText('Unable to parse file for uniprot IDs. Please reattempts with the appropriate header column labelled as uniprot_id,\nor without a header as a .txt file.\nCannot run SAFE until proper file format uploaded.')
            self.init_feat.setEnabled(False)
            df = None
            self.uids = None
        
        self.safe_run()
        
        
    def safe_run(self):
        ## -- Map uniprot IDs to GO Annotations (make matrix) -- ##
        # Use the GO to uniprot ID file already made
        goa = pd.read_csv('./runfiles/go_annotations.csv', index_col=0)

        # Now read in our file with our uniprot IDs
        uids = pd.DataFrame(self.uids, columns=['ACC_ID'])

        # remove -# from some UIDs
        uids.loc[uids['ACC_ID'].str.contains('-'),'ACC_ID'] = uids[uids['ACC_ID'].str.contains('-')]['ACC_ID'].str[:-2]

        # generate uid list
        uids_to_map = pd.DataFrame(uids['ACC_ID'])
        uids_to_map = uids_to_map.drop_duplicates().reset_index(drop=True)
        uids_to_map = uids_to_map.rename(columns={'ACC_ID':'uid'})

        # Now isolate uids_to_map within the goa matrix
        mapped = pd.merge(uids_to_map, goa, on='uid')

        # Now we need to make a matrix of uid x go term
        self.go_matrix = pd.crosstab(mapped['uid'], mapped['go_term'])
        
        
        ## -- load network -- ##
        sf = safe.SAFE() # create safe object
        sf.load_network(network_file = './runfiles/safe_huri_network.gpickle')
        
        ## -- load attributes -- ##
        sf.load_attributes(attribute_file = self.go_matrix)
        
        ## -- define_neighborhoods -- ##
        sf.define_neighborhoods(node_distance_metric = self.node_distance_metric, neighborhood_radius = self.neighborhood_radius)
        
        ## -- compute_pvalues -- ##
        if self.full_force:
            nprocesses = mp.cpu_count() - 1
        else:
            nprocesses = 1
        sf.compute_pvalues(multiple_testing = self.multiple_testing, background = self.background, num_processes = nprocesses)
        
        ## -- define_top_attributes -- ##
        sf.define_top_attributes()
        
        ## -- define_domains -- ##
        sf.define_domains(attribute_distance_threshold = self.attribute_distance_threshold)
        
        ## -- trim_domains -- ##
        sf.trim_domains()
        
        ## -- plot_composite_network (modified, displays legend) + save network -- ##
        filesavetime = time.strftime("%Y-%m-%d_%H.%M.%S")
        sf.plot_composite_network_with_legend(background_color = self.background_color, save_fig = self.save_fig + 'ProteoML_Genius_'+ filesavetime +'.png')
        # show user network
        self.sf_disp_filled = QPixmap(self.save_fig + 'ProteoML_Genius_'+ filesavetime +'.png')
        self.sf_disp_filled_sm = self.sf_disp_filled.scaledToWidth(800)  # scale as image is large!
        self.sf_label.setPixmap(self.sf_disp_filled_sm)
        
        ## -- show user domains in legend (GO-annotations and # of associated UIDs) -- ##
        # build df to show
        go_desc = []
        go_term = []
        go_net = []
        go_data = []
        go_solo = goa.drop_duplicates(subset='go_term').reset_index(drop=True)
        # annotate domains first with go_description
        for r, h in sf.domains.iterrows():
            if r > 0:
                labels = h['label'].split(', ')[1:]
                for i in labels:
                    # cleanup GO term
                    go_search = 'GO:' + i
                    go_term.append(go_search)
                    # find detailed description
                    go_desc.append(go_solo.loc[go_solo['go_term']==go_search]['go_description'].values[0])
                    # find # of uids associated within initial dataset
                    go_data.append(self.go_matrix[go_search].sum())
                    # find # of nodes associated mapped to network
                    go_net.append(len(sf.node2domain[sf.node2domain['primary_domain'] == h['id']]))
         
        self.mapped_goa = pd.DataFrame(zip(go_data, go_net, go_term, go_desc), columns=['Proteins in Dataset', 'Proteins in Network', 'GO Term',  'Detailed Description'])
         
        # show df to user
        self.go_model = PandasModel(self.mapped_goa)
        self.view.setModel(self.go_model)
        
        
if __name__ == '__main__':
    app = QApplication([])
    window = Window()
    window.show()
    sys.exit(app.exec())

