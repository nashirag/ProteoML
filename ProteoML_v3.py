## ProteoML Genius ##

# IMPORT RELEVANT PACKAGES
import sys
from PyQt6.QtCore import *
from PyQt6.QtWidgets import *
from PyQt6.QtGui import *
import os
import time
import shutil
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
from sklearn.exceptions import ConvergenceWarning, FitFailedWarning
# FOR FEATURE GENERATION
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import rdmolfiles
import itertools
# FOR ML
from collections import Counter
from sklearn.model_selection import train_test_split
from imblearn.pipeline import Pipeline
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score
from sklearn.metrics import auc, precision_score, recall_score, confusion_matrix, roc_curve, roc_auc_score
from numpy import mean
from sklearn.metrics import recall_score
from sklearn.metrics import make_scorer
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
# ml models
from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import ComplementNB
from sklearn import tree
from sklearn.neighbors import KNeighborsClassifier
from sklearn import svm
from sklearn.ensemble import BaggingClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import GradientBoostingClassifier
# sampling methods
from imblearn.over_sampling import RandomOverSampler, SMOTE, BorderlineSMOTE, SVMSMOTE, KMeansSMOTE, ADASYN
from imblearn.under_sampling import RandomUnderSampler, CondensedNearestNeighbour
from imblearn.under_sampling import TomekLinks, EditedNearestNeighbours
from imblearn.under_sampling import NeighbourhoodCleaningRule, OneSidedSelection
from imblearn.combine import SMOTEENN, SMOTETomek
from scipy.stats import loguniform
# FOR PYINSTALLER JOB HANDLING
import multiprocessing
import matplotlib.backends.backend_pdf

multiprocessing.freeze_support()

print('Imported libraries successfully.')

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
        
        
class Window(QScrollArea):
    def __init__(self):
        
        warnings.filterwarnings("ignore", category=FitFailedWarning)
        
        super().__init__(parent=None)
        
        wid = QWidget(self)
        self.setWidget(wid)
        self.setWidgetResizable(True)
        main_layout = QGridLayout()
        wid.setLayout(main_layout)
        
        self.setWindowTitle("ProteoML Genius v1.3")   # label window
        
        self.setMinimumSize(930, 800)
        
        #app.aboutToQuit.connect(self.closeEvent)
        
        # CREATE TAB WIDGET
        self.tab = QTabWidget()
        
        self.meta_model_exists = False
        self.base_model = None
        
        ############################################ DOWNLOAD PEPTIDE LIST TAB ############################################
        
        # Set up the tab
        peplist = QWidget(self)
        layout = QGridLayout()
        peplist.setLayout(layout)
        
        # Label telling user what's going on
        layout.addWidget(QLabel("Please select one or more of the PTM datasets below, then click download:"), 0, 0, 1, 3)
        layout.addWidget(QLabel(" "), 1, 0, 1, 3)   # spacer
        
        self.meth = QPushButton("Lysine Methylation")
        self.meth.setCheckable(True)   # allows an on/off state
        self.acet = QPushButton("Acetylation")
        self.acet.setCheckable(True)
        self.sumo = QPushButton("Sumoylation")
        self.sumo.setCheckable(True)
        self.arg_meth = QPushButton("Arginine Methylation")
        self.arg_meth.setCheckable(True)
        self.ubiq = QPushButton("Ubiquitination")
        self.ubiq.setCheckable(True)
        self.y_phos = QPushButton("Phosphorylation - Tyrosine")
        self.y_phos.setCheckable(True)
        self.s_phos = QPushButton("Phosphorylation - Serine")
        self.s_phos.setCheckable(True)
        self.t_phos = QPushButton("Phosphorylation - Threonine")
        self.t_phos.setCheckable(True)
        
        self.endtext = QLabel(" ")   # label for displaying message later
        
        self.down = QPushButton("Download File(s)")   # download button
        
        # Make a group for the dataset buttons
        self.down_group = QButtonGroup(peplist)
        self.down_group.addButton(self.meth)
        self.down_group.addButton(self.acet)
        self.down_group.addButton(self.sumo)
        self.down_group.addButton(self.arg_meth)
        self.down_group.addButton(self.ubiq)
        self.down_group.addButton(self.y_phos)
        self.down_group.addButton(self.s_phos)
        self.down_group.addButton(self.t_phos)
        self.down_group.setExclusive(False)
        
        self.down.clicked.connect(self.downstate)   # Download button clicked - initiate popup download dialogue!
        
        # add file buttons to layout
        layout.addWidget(self.meth, 2, 0)
        layout.addWidget(self.acet, 2, 1)
        layout.addWidget(self.sumo, 2, 2)
        layout.addWidget(self.arg_meth, 3, 0)
        layout.addWidget(self.ubiq, 3, 2)
        layout.addWidget(self.y_phos, 4, 1)
        layout.addWidget(self.s_phos, 4, 0)
        layout.addWidget(self.t_phos, 4, 2)
        layout.addWidget(QLabel(" "), 5, 0, 1, 3)   # spacer
        layout.addWidget(self.down, 6, 0, 1, 3)   # download button
        layout.addWidget(self.endtext, 7, 0, 2, 3)   # end text (blank at first)
        
        layout.addWidget(QLabel(" "), 8, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 9, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 10, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 11, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 12, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 13, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 14, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 15, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 16, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 17, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 18, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 19, 0, 1, 3)   # spacer
        
        ############################################ LOAD PEPTIDE RESULTS TAB ############################################
        
        # Import protdcal data for later use
        print('Reading in ProtDCal data...')
        pdcal_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'runfiles/protdcal_features.csv'))
        self.protdcal = pd.read_csv(pdcal_file, index_col=0)
        print('Complete.')
        
        # Import go annotations for later use
        print('Reading in GO Annotations...')
        goa_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'runfiles/go_annotations.csv'))
        self.goa = pd.read_csv(goa_file, index_col=0)
        print('Complete.')
        
        # Set up the tab
        upload = QWidget(self)
        up_layout = QGridLayout()
        upload.setLayout(up_layout)
        
        # Label telling user what's going on
        up_layout.addWidget(QLabel("Please upload your PTM results:"), 0, 0)
        
        # Upload dataset Button
        self.up = QPushButton("Upload Result Dataset")
        up_layout.addWidget(self.up, 1, 0)
        
        self.user_feat = None
        self.up.clicked.connect(self.upstate)   # upload button pushed! Initiate dialogue.
        
        
        self.up_endtext = QLabel(" ")   # label for displaying message later
        up_layout.addWidget(self.up_endtext, 2, 0)
        
        up_layout.addWidget(QLabel(" "), 3, 0)   # spacer
        
        # Initiate feature generation Button
        self.init_feat = QPushButton("Initiate Feature Generation")
        self.init_feat.setEnabled(False)
        up_layout.addWidget(self.init_feat, 4, 0)
        
        self.feat_x = None
        self.feat_y = None
        self.mapped_goa = pd.DataFrame([['GO:0000000',0,'None']], columns=['go_term', 'go_count', 'go_description'])
        
        self.init_feat.clicked.connect(lambda: self.feature_generate())
        
        
        # Progress bar for feature generation
        self.progress = QProgressBar(self)
        up_layout.addWidget(self.progress, 5, 0)
        
        self.post_feat = QLabel(" ")   # label for displaying message later
        up_layout.addWidget(self.post_feat, 6, 0)
        
        layout.addWidget(QLabel(" "), 7, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 8, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 9, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 10, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 11, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 12, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 13, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 14, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 15, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 16, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 17, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 18, 0, 1, 3)   # spacer
        layout.addWidget(QLabel(" "), 19, 0, 1, 3)   # spacer
        
        
        
        ############################################ FIT MODELS TAB ############################################
        
        # Set up the tab
        models = QWidget(self)
        mod_layout = QGridLayout()
        models.setLayout(mod_layout)
        
        mod_layout.addWidget(QLabel("Click the Button Below to Generate a Model for the Data:"), 0, 0, 1, 2)
        
        # Button for user to sort through models and sampling methods - will return best fit
        self.autofit = QPushButton("Fit Models and Balance Data - Return Best Fit")
        mod_layout.addWidget(self.autofit, 1, 0, 1, 2)
        
        # Progress bar for auto model fitting
        self.auto_modfit_progress = QProgressBar(self)
        mod_layout.addWidget(self.auto_modfit_progress, 2, 0, 1, 2)
        
        self.autofit.clicked.connect(self.model_fitting)   # autofit button clicked - initiate ml autofitting
        
        # Drop-down menus for models and sampling methods
        mod_layout.addWidget(QLabel("Machine Learning Model:"), 3, 0, 1, 1)
        self.modselect = QComboBox()
        self.modselect.addItem("Dummy Classifier")
        self.modselect.addItem("Logistic Regression")
        self.modselect.addItem("Linear Discriminant Analysis")
        self.modselect.addItem("Complement Naive Bayes")
        self.modselect.addItem("Decision Tree Classifier")
        self.modselect.addItem("K-Nearest Neighbours Classifier (k-NN)")
        self.modselect.addItem("Support Vector Classifier (SVC)")
        self.modselect.addItem("Bagging Classifier")
        self.modselect.addItem("Random Forest Classifier")
        self.modselect.addItem("Extra Trees Classifier")
        self.modselect.addItem("Gradient Boosting Classifier")
        mod_layout.addWidget(self.modselect, 4, 0)
        
                
        mod_layout.addWidget(QLabel("Data Balancing Method:"), 3, 1, 1, 1)
        self.sampselect = QComboBox()
        self.sampselect.addItem("None")
        self.sampselect.addItem("Oversampling: Random Oversampler")
        self.sampselect.addItem("Oversampling: SMOTE")
        self.sampselect.addItem("Oversampling: Borderline SMOTE")
        self.sampselect.addItem("Oversampling: SVMSMOTE")
        #self.sampselect.addItem("Oversampling: K Means SMOTE")
        self.sampselect.addItem("Oversampling: ADASYN")
        self.sampselect.addItem("Undersampling: Random Undersampler")
        self.sampselect.addItem("Undersampling: Condensed Nearest Neighbour")
        self.sampselect.addItem("Undersampling: Tomek Links")
        self.sampselect.addItem("Undersampling: Edited Nearest Neighbour")
        self.sampselect.addItem("Undersampling: Neighbourhood Cleaning Rule")
        self.sampselect.addItem("Undersampling: One Sided Selection")
        self.sampselect.addItem("Combined: SMOTE-ENN")
        self.sampselect.addItem("Combined: SMOTE Tomek")
        mod_layout.addWidget(self.sampselect, 4, 1)
        
        self.f1_score_label = QLabel("Current F1 Score: ")   # label for displaying f1 score later
        mod_layout.addWidget(self.f1_score_label, 5, 0, 1, 2)
        
        # Button for user to regenerate model fitting with their selection
        self.userfit = QPushButton("Re-Run Fitting with New Selection")
        mod_layout.addWidget(self.userfit, 6, 0, 1, 2)
        
        self.userfit.clicked.connect(self.model_fitting)   # autofit button clicked - initiate ml autofitting
        
        # Progress bar for auto model fitting
        self.user_modfit_progress = QProgressBar(self)
        mod_layout.addWidget(self.user_modfit_progress, 7, 0, 1, 2)
        
        #mod_layout.addWidget(QLabel(" "), 7, 0, 1, 2)   # spacer
        
        # Display for pdfs of PR and ROC curves
        mod_layout.addWidget(QLabel("Precision-Recall Curve"), 8, 0, 1, 2)
        pr_disp_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "runfiles/pr_blank.png"))
        self.pr_disp = QPixmap(pr_disp_file)
        self.pr_label = QLabel()
        self.pr_label.setPixmap(self.pr_disp)
        mod_layout.addWidget(self.pr_label, 9, 0)

        mod_layout.addWidget(QLabel("ROC Curve"), 8, 1, 1 , 2)
        roc_disp_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "runfiles/roc_blank.png"))
        self.roc_disp = QPixmap(roc_disp_file)
        self.roc_label = QLabel()
        self.roc_label.setPixmap(self.roc_disp)
        mod_layout.addWidget(self.roc_label, 9, 1)
        
        
        # Labels for metrics below PR and ROC curves
        mod_layout.addWidget(QLabel(" - Metrics at Specified Threshold - "), 10, 0, 1, 1)
        #mod_layout.addWidget(QLabel("Threshold: "), 11, 0, 1, 1)
        self.user_thresh = QLineEdit('0.5')
        mod_layout.addWidget(self.user_thresh, 11, 0, 1, 1)
        self.sens_lab = QLabel("Sensitivity:  ")
        self.spec_lab = QLabel("Specificity: ")
        self.prec_lab = QLabel("Precision: ")
        self.rec_lab = QLabel("Recall: ")
        mod_layout.addWidget(self.prec_lab, 12, 0)
        mod_layout.addWidget(self.rec_lab, 13, 0)
        mod_layout.addWidget(self.sens_lab, 12, 1)
        mod_layout.addWidget(self.spec_lab, 13, 1)
        
        mod_layout.addWidget(QLabel(""), 14, 0, 1, 2)
        
        # Space for commonly associated GO annotations
        self.go = QLabel(" - Highly Associated GO Terms - ")
        mod_layout.addWidget(self.go, 15, 0, 1, 1)
        self.go_model = PandasModel(self.mapped_goa)
        self.view = QTableView()
        self.view.setModel(self.go_model)
        self.view.horizontalHeader().resizeSection(0, 150)
        self.view.horizontalHeader().resizeSection(1, 150)
        self.view.horizontalHeader().resizeSection(2, 500)
        self.view.setMinimumSize(800, 300)
        mod_layout.addWidget(self.view, 16, 0, 10, 2)
        
        self.fitwidth = models.frameGeometry().width()
        
        ############################################ METALEARNING TAB ############################################
        
        # Set up the tab
        meta = QWidget(self)
        meta_layout = QGridLayout()
        meta.setLayout(meta_layout)
        
        meta_layout.addWidget(QLabel("This step is optional. \nTo create a MetaLearning Model, the SECONDARY_ML_SCORE column in your feature dataset will be used with\nthose of your base predictor (created on the previous tab)."), 0, 0, 1, 2)
        
        self.meta_secondary = None
        
        # Output for meta learning model dataset
        self.meta_data_confirm = QLabel(" ")   # Output after features have been read in
        meta_layout.addWidget(self.meta_data_confirm, 1, 0, 1, 2)
        
        # Button for user to sort through models and sampling methods - will return best fit
        self.meta_autofit = QPushButton("Fit Models and Balance Data - Return Best Fit")
        meta_layout.addWidget(self.meta_autofit, 3, 0, 1, 2)
        
        # Progress bar for auto model fitting
        self.auto_metafit_progress = QProgressBar(self)
        meta_layout.addWidget(self.auto_metafit_progress, 4, 0, 1, 2)
        
        self.meta_autofit.clicked.connect(self.meta_model_fitting)   # autofit button clicked - initiate ml autofitting
        
        # Drop-down menus for models and sampling methods
        meta_layout.addWidget(QLabel("Machine Learning Model:"), 5, 0, 1, 1)
        self.meta_modselect = QComboBox()
        self.meta_modselect.addItem("Dummy Classifier")
        self.meta_modselect.addItem("Logistic Regression")
        self.meta_modselect.addItem("Linear Discriminant Analysis")
        self.meta_modselect.addItem("Complement Naive Bayes")
        self.meta_modselect.addItem("Decision Tree Classifier")
        self.meta_modselect.addItem("K-Nearest Neighbours Classifier (k-NN)")
        self.meta_modselect.addItem("Support Vector Classifier (SVC)")
        self.meta_modselect.addItem("Bagging Classifier")
        self.meta_modselect.addItem("Random Forest Classifier")
        self.meta_modselect.addItem("Extra Trees Classifier")
        self.meta_modselect.addItem("Gradient Boosting Classifier")
        meta_layout.addWidget(self.meta_modselect, 6, 0)
        
                
        meta_layout.addWidget(QLabel("Data Balancing Method:"), 5, 1, 1, 1)
        self.meta_sampselect = QComboBox()
        self.meta_sampselect.addItem("None")
        self.meta_sampselect.addItem("Oversampling: Random Oversampler")
        self.meta_sampselect.addItem("Oversampling: SMOTE")
        self.meta_sampselect.addItem("Oversampling: Borderline SMOTE")
        self.meta_sampselect.addItem("Oversampling: SVMSMOTE")
        #self.meta_sampselect.addItem("Oversampling: K Means SMOTE")
        self.meta_sampselect.addItem("Oversampling: ADASYN")
        self.meta_sampselect.addItem("Undersampling: Random Undersampler")
        self.meta_sampselect.addItem("Undersampling: Condensed Nearest Neighbour")
        self.meta_sampselect.addItem("Undersampling: Tomek Links")
        self.meta_sampselect.addItem("Undersampling: Edited Nearest Neighbour")
        self.meta_sampselect.addItem("Undersampling: Neighbourhood Cleaning Rule")
        self.meta_sampselect.addItem("Undersampling: One Sided Selection")
        self.meta_sampselect.addItem("Combined: SMOTE-ENN")
        self.meta_sampselect.addItem("Combined: SMOTE Tomek")
        meta_layout.addWidget(self.meta_sampselect, 6, 1)
        
        self.meta_f1_score_label = QLabel("Current F1 Score: ")   # label for displaying f1 score later
        meta_layout.addWidget(self.meta_f1_score_label, 7, 0, 1, 2)
        
        # Button for user to regenerate model fitting with their selection
        self.meta_userfit = QPushButton("Re-Run Fitting with New Selection")
        meta_layout.addWidget(self.meta_userfit, 8, 0, 1, 2)
        
        self.meta_userfit.clicked.connect(self.meta_model_fitting)   # autofit button clicked - initiate ml autofitting
        
        # Progress bar for auto model fitting
        self.meta_user_modfit_progress = QProgressBar(self)
        meta_layout.addWidget(self.meta_user_modfit_progress, 9, 0, 1, 2)
        
        #mod_layout.addWidget(QLabel(" "), 7, 0, 1, 2)   # spacer
        
        # Display for pdfs of PR and ROC curves
        meta_layout.addWidget(QLabel("Precision-Recall Curve"), 10, 0, 1, 2)
        self.meta_pr_disp = QPixmap(pr_disp_file)
        self.meta_pr_label = QLabel()
        self.meta_pr_label.setPixmap(self.meta_pr_disp)
        meta_layout.addWidget(self.meta_pr_label, 11, 0)

        meta_layout.addWidget(QLabel("ROC Curve"), 10, 1, 1 , 2)
        self.meta_roc_disp = QPixmap(roc_disp_file)
        self.meta_roc_label = QLabel()
        self.meta_roc_label.setPixmap(self.meta_roc_disp)
        meta_layout.addWidget(self.meta_roc_label, 11, 1)
        
        # Labels for metrics below PR and ROC curves
        meta_layout.addWidget(QLabel(" - Metrics at a Threshold of 0.5 - "), 12, 0, 1, 1)
        self.meta_sens_lab = QLabel("Sensitivity:  ")
        self.meta_spec_lab = QLabel("Specificity: ")
        self.meta_prec_lab = QLabel("Precision: ")
        self.meta_rec_lab = QLabel("Recall: ")
        meta_layout.addWidget(self.meta_prec_lab, 13, 0)
        meta_layout.addWidget(self.meta_rec_lab, 14, 0)
        meta_layout.addWidget(self.meta_sens_lab, 13, 1)
        meta_layout.addWidget(self.meta_spec_lab, 14, 1)
        
        ############################################ EXPERIMENTAL PREDICTIONS TAB ############################################
        # 1A - Select file + upload
                
        # Set up the tab
        experimental = QWidget(self)
        exp_layout = QGridLayout()
        experimental.setLayout(exp_layout)
        
        # Label telling user what's going on
        exp_layout.addWidget(QLabel("Now that the models have been initialized, you can insert experimental sequences for Machine Learning predictions!"), 0, 0, 1, 2)

        self.exp_set = None
        self.user_exp_set = None
        
        # 1 - Select from pregenerated feature set
        # drop down
        exp_layout.addWidget(QLabel("Option 1 - Select a Dataset:"), 1, 0, 1, 2)
        self.exp_dataset = QComboBox()
        self.exp_dataset.addItem("Surface Exposed Sites - Lysine")   #0
        self.exp_dataset.addItem("Surface Exposed Sites - Arginine")   #1
        self.exp_dataset.addItem("Surface Exposed Sites - Tyrosine")   #2
        self.exp_dataset.addItem("Surface Exposed Sites - Threonine")   #3
        exp_layout.addWidget(self.exp_dataset, 2, 0, 1, 2)
        
        # button to confirm
        self.exp_confirm = QPushButton("Download Selected Dataset (alter and/or insert scores for Meta Learning)")
        exp_layout.addWidget(self.exp_confirm, 3, 0, 1, 1)
        self.exp_proceed = QPushButton("Proceed with Selected Dataset (Basic Model only)")
        exp_layout.addWidget(self.exp_proceed, 3, 1, 1, 1)
        
        self.exp_confirm.clicked.connect(self.downstate)   # want to download file
        self.exp_proceed.clicked.connect(self.exp_proceed_upload) # want to proceed with selection - initiate feature generation
        
        exp_layout.addWidget(QLabel(" "), 4, 0)   # spacer
        
        # Upload dataset Button
        exp_layout.addWidget(QLabel("Option 2 - Upload a Dataset:"), 5, 0, 1, 2)
        self.exp_upload = QPushButton("Upload or Reupload a Dataset")
        exp_layout.addWidget(self.exp_upload, 6, 0, 1, 2)
        
        self.exp_upload.clicked.connect(self.exp_upstate)   # upload button pushed! Initiate dialogue.
        
        self.user_experimental_label = QLabel(" ")   # label for displaying message later
        exp_layout.addWidget(self.user_experimental_label, 7, 0, 1, 2)
        
        exp_layout.addWidget(QLabel(" "), 8, 0)   # spacer
        
        # 2 - Generate features
        
        # Label for showing feature generation info later
        self.feature_generating_label = QLabel(" ")
        exp_layout.addWidget(self.feature_generating_label, 9, 0)
        
        # Progress bar for feature generation
        self.exp_feat_progress = QProgressBar(self)
        exp_layout.addWidget(self.exp_feat_progress, 10, 0, 1, 2)
        
        self.exp_post_feat = QLabel(" ")   # label for displaying message later
        exp_layout.addWidget(self.exp_post_feat, 11, 0, 1, 2)
        
        exp_layout.addWidget(QLabel(" "), 12, 0, 1, 2)   # spacer
        
        # Select ML model radiobutton - select only one at a time
        self.model_group = QButtonGroup()
        self.b1 = QRadioButton("Basic Machine Learning Model")
        self.model_group.addButton(self.b1)
        self.b1.setChecked(True)
        self.b2 = QRadioButton("Meta Machine Learning Model")
        self.model_group.addButton(self.b2)
        exp_layout.addWidget(self.b1, 13, 0)
        exp_layout.addWidget(self.b2, 13, 1)
        self.b2.setCheckable(False)  # disable selecting this option until user supplies dset and features are generated
        
        # Fit dataset to ML model button
        self.exp_init_run = QPushButton("Initiate Machine Learning Run")
        self.exp_init_run.setEnabled(False)   # false until feature set selected + processed
        exp_layout.addWidget(self.exp_init_run, 14, 0, 1, 2)
        
        self.exp_init_run.clicked.connect(self.experimental_ml_run)   # initiate ML run!
        
        # Progress bar for ml run
        self.exp_ml_progress = QProgressBar(self)
        exp_layout.addWidget(self.exp_ml_progress, 15, 0, 1, 2)
        
        self.exp_post_run = QLabel(" ")   # label for displaying scored message later
        exp_layout.addWidget(self.exp_post_run, 16, 0, 1, 2)
        
        # Save score file
        self.exp_save = QPushButton("Save Score File Locally")
        self.exp_save.setEnabled(False)
        exp_layout.addWidget(self.exp_save, 17, 0, 1, 2)
        
        self.exp_post_save = QLabel(" ")   # label for displaying saved file message later
        exp_layout.addWidget(self.exp_post_save, 18, 0, 1, 2)
        
        self.launch_vis = QPushButton("Launch ProteoML Visualizer")
        self.launch_vis.setEnabled(False)
        exp_layout.addWidget(self.launch_vis, 19, 0, 2, 2)
        
        self.exp_save.clicked.connect(self.experimental_ml_run_save)  # initiate dialogue for saving scores
        #self.launch_vis.clicked.connect(self.launch_visualiser)   # launch ML experimental visualiser with dataset from predictions
        
        
        ############################################ CREATE OUR TABS ############################################
        
        # create tabs themselves
        self.tab.addTab(peplist, "1 - Download Peptide List")
        self.tab.addTab(upload, "2 - Upload Peptide Results")
        self.tab.addTab(models, "3 - Fit Models")
        self.tab.addTab(meta, "4 - Meta Learning")
        self.tab.addTab(experimental, "5 - Experimental Predictions")
        
        main_layout.addWidget(self.tab, 0, 0)
    
    def closeEvent(self, event):
        reply = QMessageBox()
        #('Message',
        reply.setIcon(QMessageBox.Icon.Warning)
        reply.setStandardButtons(QMessageBox.StandardButton.Cancel | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Yes)
        reply.setInformativeText("You're about to close the session. Would you like to save your progress?")
        reply.setDetailedText("Click cancel to return to active session. \nClick Yes to save session, or click No to ignore this warning and close the active session.")
        ret = reply.exec()
        
        if ret == QMessageBox.StandardButton.Yes:
            print('User hit yes - saving data files...')
            # Initiate save dialogue
            repo_dest = QFileDialog.getExistingDirectory(self, 'Save Session Files To')
            filesavetime = time.strftime("%Y-%m-%d_%H.%M.%S")
            repo_dest_f = str(repo_dest) + "/ProteoML_Genius_session_" + filesavetime
            
            # Save Feature Sets
            if hasattr(self, 'feat_x') or self.exp_set != None:
                # Make directory to save files in - user would have to have features if they have graphs or an experimental score result file
                os.mkdir(repo_dest_f)
                to_save_features = [self.feat_x, self.exp_set]
                to_save_features_n = ['user_defined_model_fitting_features.csv', 'user_defined_experimental_set_features.csv']
                j = 0
                for i in to_save_features:
                    if len(i) != 0:
                        i.to_csv(repo_dest_f + '/' + to_save_features_n[j])
                    else:
                        print(str(i), 'was none')
                    j += 1
            else:
                print("none of the features were populated")
            
            # Save Plot .pdfs
            if hasattr(self, 'i_pr') or hasattr(self, 'i_roc') or hasattr(self, 'm_pr') or hasattr(self, 'm_roc'):
                to_save_graphs = [self.i_pr, self.i_roc, self.m_pr, self.m_roc]
                to_save_graphs_n = ['ml_pr_curve.pdf', 'ml_roc_curve.pdf', 'meta_ml_pr_curve.pdf', 'meta_ml_roc_curve.pdf']
                j = 0
                for i in to_save_graphs:
                    if len(i) != 0:
                        shutil.copyfile(i, repo_dest_f + '/' + to_save_graphs_n[j])
                    else:
                        print(str(i), 'was none')
                    j += 1
            else:
                print("none of the graph images were populated")
            
            # Save Experimental Score File
            if hasattr(self, 'experimental_scores'):
                self.experimental_scores.to_csv(repo_dest_f + '/' + 'user_defined_experimental_set_results.csv')
            else:
                print("experimental scores weren't yet populated")
                
            event.accept()  # now we accept user's wish to close application
            
        elif ret == QMessageBox.StandardButton.No:
            print('User hit no - not saving data files, proceeding with quitting...')
            event.accept()
        else:
            print('User hit cancel - not saving data files, going back to main window...')
            event.ignore()
    
    
    def downstate(self):
        destination = QFileDialog.getExistingDirectory(self, 'Save Files To')
        if self.tab.currentIndex() == 0:
            # We're on the feature dset download page
            # DOWNLOAD (MOVE) SEQUENCE DATA TO USER SPECIFIED LOCATION
            if self.down_group.checkedButton() is None:
                self.endtext.setText("No datasets selected. Please select dataset(s) to download or proceed to next stage.")
            else:
                if self.meth.isChecked():
                    meth_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'datasets/formatted_files/excel/Lysine_Methylation_sites.xlsx'))
                    shutil.copyfile(meth_file, destination+'/Lysine_Methylation_sites.xlsx')
                if self.acet.isChecked():
                    acet_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'datasets/formatted_files/excel/Acetylation_sites.xlsx'))
                    shutil.copyfile(acet_file, destination+'/Acetylation_sites.xlsx')
                if self.t_phos.isChecked():
                    t_phos_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'datasets/formatted_files/excel/T_Phosphorylation_sites.xlsx'))
                    shutil.copyfile(t_phos_file, destination+'/T_Phosphorylation_sites.xlsx')
                if self.s_phos.isChecked():
                    s_phos_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'datasets/formatted_files/excel/S_Phosphorylation_sites.xlsx'))
                    shutil.copyfile(s_phos_file, destination+'/S_Phosphorylation_sites.xlsx')
                if self.y_phos.isChecked():
                    y_phos_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'datasets/formatted_files/excel/Y_Phosphorylation_sites.xlsx'))
                    shutil.copyfile(y_phos_file, destination+'/Y_Phosphorylation_sites.xlsx')
                if self.sumo.isChecked():
                    sumo_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'datasets/formatted_files/excel/Sumoylation_sites.xlsx'))
                    shutil.copyfile(sumo_file, destination+'/Sumoylation_sites.xlsx')
                if self.arg_meth.isChecked():
                    arg_meth_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'datasets/formatted_files/excel/Arginine_Methylation_sites.xlsx'))
                    shutil.copyfile(arg_meth_file, destination+'/Arginine_Methylation_sites.xlsx')
                if self.ubiq.isChecked():
                    ubiq_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'datasets/formatted_files/excel/Ubiquitination_sites.xlsx'))
                    shutil.copyfile(ubiq_file, destination+'/Ubiquitination_sites.xlsx')
                # confirm with user where the file ended up
                self.endtext.setText('Downloaded files to ' + destination + ' folder successfully.')
        else:
            # GENERATE AND IMPLEMENT THESE DATASETS!!
            # We're on the experimental dset download page
            desired_exp = self.exp_dataset.currentIndex()
            if desired_exp == 0:
            # Lysines
                print("User selected the surface exposed lysine dset")
            elif desired_exp == 1:
            # Arginine
                print("User selected the surface exposed arginine dset")
            elif desired_exp == 2:
            # Tyrosine
                print("User selected the surface exposed tyrosine dset")
            elif desired_exp == 3:
            # Threonine
                print("User selected the surface exposed threonine dset")
        
    def upstate(self):
        # UPLOAD FEATURES FROM USER
        upload_from = QFileDialog.getOpenFileName(self, "Upload PTM Experimental Data File")
        upload_from_dir = upload_from[0]
        
        if ".xls" in upload_from_dir:
            # Excel file
            df = pd.read_excel(upload_from_dir)
            self.up_endtext.setText('Uploaded excel file: ' + upload_from_dir)
            self.init_feat.setEnabled(True)
        elif ".csv" in upload_from_dir:
            # csv file
            df = pd.read_csv(upload_from_dir, index_col=0)
            self.up_endtext.setText('Uploaded csv file: ' + upload_from_dir)
            self.init_feat.setEnabled(True)
        elif ".tsv" in upload_from_dir:
            df = pd.read_csv(upload_from_dir, sep='\t', index_col=0)
            self.up_endtext.setText('Uploaded tsv file: ' + upload_from_dir)
            self.init_feat.setEnabled(True)
        else:
            # improper format >:(
            self.up_endtext.setText('Unable to upload file. Please attempt with a .xlsx, .csv, or .tsv file.')
            self.init_feat.setEnabled(False)
            df = 'NaN'
        self.user_feat = df

    def feature_generate(self):
        # DOWNSTATE FOR FEATURE GENERATE BUTTON - INITIATE FEATURE GENERATION SCRIPT FOR PASSED FILE
        
        # Format user defined features for FeatureGen call
        temp_annotation = self.user_feat['MOD_RSD'].str.split("-", expand=True)[0]
        temp_annotation = temp_annotation.str[1:]
        uid_annotation = self.user_feat['ACC_ID'] + "_" + temp_annotation
        self.user_feat = self.user_feat.set_index(uid_annotation, drop=True)
        
        # Remove duplicates
        feat_dup_len = (len(self.user_feat))
        self.user_feat = self.user_feat[~self.user_feat.index.duplicated(keep='first')]
        feat_dup_len = feat_dup_len - (len(self.user_feat))
        
        # Pull sequences
        sequences = self.user_feat['SITE_+/-7_AA']
        
        # Pull user defined features' uniprot IDs from GO table - generate associated GO annotations
        feat_uids = pd.DataFrame(self.user_feat['ACC_ID'])
        feat_uids = feat_uids.rename(columns={'ACC_ID':'uid'})
        feat_uids = feat_uids.drop_duplicates().reset_index(drop=True)
        print(feat_uids.columns)
        goa_uids = pd.merge(feat_uids, self.goa, on='uid')
        goa_counts = pd.DataFrame(goa_uids['go_term'].value_counts())
        goa_counts = goa_counts.reset_index(drop=False).rename(columns={'index':'go_term', 'count':'go_count'})
        goa_unique = goa_uids.drop_duplicates(subset='go_term')
        self.mapped_goa = pd.merge(goa_counts, goa_unique, on='go_term').drop(columns=['uid', 'gene', 'go_association', 'go_ref', 'go_b_p_f'])   # finalised df with counts to be displayed to user
        
        # Update dataframe on model fitting tab
        self.go_model = PandasModel(self.mapped_goa)
        self.view.setModel(self.go_model)
        
        # Now get into generating our features!
        
        # Create df for results to go into
        v, h = FeatureGen(sequences[0], self.protdcal)
        features = pd.DataFrame(columns=h)
        features.loc[len(features)] = v
        
        increment_by = 100/len(sequences)
        completed = increment_by
        i = 1
        
        # Go through rest of sequences to generate feature set
        while i < len(sequences):
            ts = sequences[i]
            value, header = FeatureGen(ts, self.protdcal)
            features.loc[len(features)] = value
            completed+=increment_by
            i+=1
            QApplication.processEvents()
            self.progress.setValue(completed)
            self.post_feat.setText("Feature Generation is " + str(round(completed, 1)) + "% complete.")

        # Make the index the same as our initial dataframe
        self.feat_x = features.set_index(self.user_feat.index)
        
        # Isolate the methylated condition from the sequences as our y value
        self.feat_y = self.user_feat['EXPERIMENTALLY_ACTIVE']
        
        self.progress.setValue(100)
        
        if 'SECONDARY_ML_SCORE' not in self.user_feat.columns or (self.user_feat['SECONDARY_ML_SCORE'] == 0).all():
            self.meta_data_confirm.setText("Cannot proceed with Meta Learning as the SECONDARY_ML_SCORE column in the feature dataset doesn't contain any scores.\nTo proceed, populate the SECONDARY_ML_SCORE column with numeric scores from another computational method.")
        else:
            self.meta_data_confirm.setText("SECONDARY_ML_SCORE successfully uploaded. Proceed with Meta Learning!")
            self.meta_secondary = self.user_feat['SECONDARY_ML_SCORE']
            self.ms_full = self.user_feat[['SECONDARY_ML_SCORE', 'ACC_ID', 'SITE_LOC']]
        
        self.post_feat.setText("Feature Generation has finished! \nThe chosen dataset has " + str(Counter(self.feat_y)[1]) + " positives and " + str(Counter(self.feat_y)[0]) + " negatives.\n"+str(feat_dup_len)+" duplicated features removed.\nYou may proceed with Model Fitting by clicking the next tab.")
        
        
    def model_fitting(self):
        # Ignore warnings - need to find a better way of suppressing output later on
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        
        # Get button pushed on tab - tells us which ML fitting mode we're entering
        sender = self.sender()
        
        # CREATE MODELS + SAMPLING METHODS FOR LATER USE
        # create our models
        dummy = DummyClassifier(strategy='most_frequent')
        logit = LogisticRegression(max_iter=100)
        lda = LinearDiscriminantAnalysis()
        nb = ComplementNB()
        dtree = tree.DecisionTreeClassifier()
        knn = KNeighborsClassifier()
        svc = svm.SVC()
        bagged = BaggingClassifier()
        rand_forest = RandomForestClassifier()
        ext_trees = ExtraTreesClassifier()
        gboost = GradientBoostingClassifier()
        
        # create list of our models to loop through
        models = [dummy, logit, lda, nb, dtree, knn, svc, bagged, rand_forest, ext_trees, gboost]
        
        # create our sampling methods
        # oversamplers
        over = RandomOverSampler()
        smote = SMOTE()
        border_smote = BorderlineSMOTE()
        svm_smote = SVMSMOTE()
        km_smote = KMeansSMOTE()
        adasyn = ADASYN()
        # undersamplers
        under = RandomUnderSampler()
        cnn = CondensedNearestNeighbour()
        tomek = TomekLinks()
        enn = EditedNearestNeighbours()
        n_cleaning = NeighbourhoodCleaningRule()
        onesided = OneSidedSelection()
        # combined samplers
        smoteenn = SMOTEENN()
        smotetomek = SMOTETomek()
        # none
        none = "none"
        
        # create list of sampling methods to loop through - removed ->
        samplings = [none, over, smote, border_smote, svm_smote, km_smote,
                     adasyn, under, cnn, tomek, enn, n_cleaning,
                     onesided, smoteenn, smotetomek]
                     
                     
        # NOW ONTO MODEL FITTING
        
        # AUTO-FITTING MODEL + SAMPLING METHOD
        if sender.text() == 'Fit Models and Balance Data - Return Best Fit':
        
            print('Onto the auto-model fitting.')   # db
            
            ## FIT VARIOUS MODELS AND SAMPLING METHODS - SELECT FOR THE BEST ONE :) ##
            
            f1 = []   # scoring metric
            increment_by = (100/3)/len(models)   # progress bar
            completed = increment_by   # update progress bar
            
            # First loop through models list and generate f1 scoring
            for m in models:
                print('Basic model fitting')   #db
                # We're on the basic model - use the simple metric
                metric = k_fold(self.feat_x, self.feat_y, m, 'none', 'none')
                completed+=increment_by
                QApplication.processEvents()
                self.auto_modfit_progress.setValue(round(completed))
                f1.append(metric)

            
            # Return model with biggest f1 score
            model_max_f1 = max(f1)
            model_max = models[f1.index(max(f1))]

            
            # Secondly we go through sampling methods
            
            f1_score = 0
            f1_s = []
            increment_by = (100/3)/len(samplings)
            print('now onto sampling methods')   #db
            
            for s in samplings:
                if s == none:
                    metric = model_max_f1
                else:
                    print('Basic model fitting')   #db
                    # We're on the basic model - use the simple metric
                    metric = k_fold(self.feat_x, self.feat_y, m, s, 'sampling')
                    completed+=increment_by
                    QApplication.processEvents()
                    self.auto_modfit_progress.setValue(round(completed))
                f1_s.append(metric)

            
            # Identify sampling method with biggest f1 score
            samp_max_f1 = max(f1_s)
            samp_max = samplings[f1_s.index(max(f1_s))]
            
            # Oversample x and y datasets and define best performing model
            if samp_max == none:
                sampler = None
            else:
                sampler = samp_max
            final_model = model_max

            # Return final f1 score
            f1_score = samp_max_f1
            
            # Return strings of model and sampling
            model_str = str(model_max)
            sampling_str = str(samp_max)
            
            model_index = f1.index(max(f1))
            sampling_index = f1_s.index(max(f1_s))
            
            print('best model', model_str)   #db
            print('best sampling method', sampling_str)   #db
                        
            # Update drop down menus with selected model and sampling methods
            self.modselect.setCurrentIndex(model_index)
            self.sampselect.setCurrentIndex(sampling_index)

        
        # USER-SELECTED MODEL + SAMPLING METHOD
  
        elif sender.text() == 'Re-Run Fitting with New Selection':
            
            # USER SELECTED THEIR OWN MODEL + SAMPLING COMBINATION
            QApplication.processEvents()   # update progress bar
            self.user_modfit_progress.setValue(15)
            
            # user selected model
            i_u_mod = self.modselect.currentIndex()
            # user selected sampling method
            i_u_samp = self.sampselect.currentIndex()
            
            # user-defined pipeline
            if samplings[i_u_samp] == none:
                sampler = None
                final_model = models[i_u_mod]
                f1_score = k_fold(self.feat_x, self.feat_y, final_model, 'none', 'none')
            else:
                final_model = models[i_u_mod]
                sampler = samplings[i_u_samp]
                f1_score = k_fold(self.feat_x, self.feat_y, final_model, sampler, 'sampling')
            
            # Return strings of model and sampling
            model_str = str(models[i_u_mod])
            sampling_str = str(samplings[i_u_samp])
                
        
        # HYPERPARAMETER SEARCH FOR MODEL OF BEST FIT
        if sampler != None:
            # Apply sampling method to training data
            sampled_x, sampled_y = sampler.fit_resample(self.feat_x, self.feat_y)
        else:
            sampled_x = self.feat_x
            sampled_y = self.feat_y
        
        if self.modselect.currentIndex() != 0:
            # Be sure we aren't using the dummy classifier (no hyperparameters for that)
            params = tuning(self.modselect.currentIndex(), final_model, sampled_x, sampled_y)
        
            # Apply hyperparameters
            final_model.set_params(**params)
            print(final_model.get_params())
        
        self.base_model = final_model
        
        ## GENERATE PR AND ROC CURVE
        
        # FOR THE BASIC MODEL
        
        threshold = 0.0
        rounds = 3
        count = 0
        t_rec = []
        rec = []
        t_prec = []
        prec = []
        t_spec = []
        spec = []
        t_sens = []
        sens = []
        thresh = []
        
        increment_by = (100/3)/(100*rounds)
        user_increment_by = 70/(100*rounds)
        user_completed = user_increment_by + 30
        
        
        while count < rounds:
            while threshold <= 1:

                # split into test + train sets
                X_train, X_test, y_train, y_test = train_test_split(self.feat_x
                                                                    , self.feat_y
                                                                    , stratify=self.feat_y)
                
                if sampler != None:
                    # Apply sampling method to training data
                    X_train, y_train = sampler.fit_resample(X_train, y_train)
                
                # fit the model to training data
                final_model.fit(X_train, y_train)

                # calculate probability on testing data
                y_proba = final_model.predict_proba(X_test)
                y_proba = y_proba[:, [1]]   #select the probability for the positive case only

                # select just the scores over our threshold
                cond = y_proba >= threshold
                # convert to 1 or 0 values
                y_pred = np.where((y_proba>=threshold), 1, y_proba)
                y_pred = np.where((y_pred<threshold), 0, y_pred)
                # calculate recall
                recall = recall_score(y_test, y_pred)
                t_rec.append(recall)
                # calculate precision
                precision = precision_score(y_test, y_pred)
                t_prec.append(precision)
                # calculate specificity
                tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
                specificity = tn/(tn+fp)
                t_spec.append(specificity)
                sensitivity = tp/(tp+fn)
                t_sens.append(sensitivity)

                # update threshold
                thresh.append(threshold)
                threshold += 0.01
                
                                    
                QApplication.processEvents()
                
                if sender.text() == 'Fit Models and Balance Data - Return Best Fit':
                    completed+=increment_by
                    self.auto_modfit_progress.setValue(round(completed))
                elif sender.text() == 'Re-Run Fitting with New Selection':
                    user_completed+=user_increment_by
                    self.user_modfit_progress.setValue(round(user_completed))
                    
            # Add metric array to metric array
            if count == 0:
                rec = t_rec
                prec = t_prec
                spec = t_spec
                sens = t_sens
            else:
                rec = [r + tr for r, tr in zip(rec, t_rec)]
                prec = [p + tp for p, tp in zip(prec, t_prec)]
                spec = [s + ts for s, ts in zip(spec, t_spec)]
                sens = [n + tn for n, tn in zip(sens, t_sens)]
            
            count+=1
            threshold=0.0
        
        rec = [r / rounds for r in rec]
        prec = [p / rounds for p in prec]
        spec = [s / rounds for s in spec]
        sens = [n / rounds for n in sens]
        
        fpr, tpr, _ = roc_curve(y_test,  y_proba)
        fpr = fpr.tolist()
        tpr = tpr.tolist()
        auc_m = roc_auc_score(y_test, y_proba)
        
        pr_m = pd.DataFrame(list(zip(thresh, rec, spec, prec, sens)), columns=['Threshold'
                                                                                  , 'Recall'
                                                                                  , 'Specificity'
                                                                                  , 'Precision'
                                                                                  , 'Sensitivity'])
        roc_m = pd.DataFrame(list(zip(fpr, tpr)), columns = ['fpr', 'tpr'])
        
        # UPDATE OUTPUT + GENERATE PLOTS
            
        # Update f1 score to user
        self.f1_score_label.setText("Current F1 Score: " + str(round(f1_score, 2)))
        
        
        # Call PR and ROC curve generation
        ms = ''.join(filter(str.isalnum, model_str))
        ss =''.join(filter(str.isalnum, sampling_str))
        save_annotation = 'runfiles/temp/Basic_ML_Model_' + ms + '_' + ss
        Plot(pr_m, (save_annotation+'_pr_curve'), roc_m, auc_m, (save_annotation+'_roc_curve'))
        
        # Update PR and ROC curves to user
        pr_png_file = os.path.abspath(os.path.join(os.path.dirname(__file__), (save_annotation+'_pr_curve.png')))
        roc_png_file = os.path.abspath(os.path.join(os.path.dirname(__file__), (save_annotation+'_roc_curve.png')))
        self.pr_label.setPixmap(QPixmap(pr_png_file))
        self.roc_label.setPixmap(QPixmap(roc_png_file))
        
        # Paths to images for file saving later
        self.i_pr = save_annotation+'_pr_curve.pdf'
        self.i_roc = save_annotation+'_roc_curve.pdf'
        
        # Update metric scores
        self.sens_lab.setText("Sensitivity: " + str(round(pr_m.iloc[50]['Sensitivity'], 2)))
        self.spec_lab.setText("Specificity: " + str(round(pr_m.iloc[50]['Specificity'], 2)))
        self.prec_lab.setText("Precision: " + str(round(pr_m.iloc[50]['Precision'], 2)))
        self.rec_lab.setText("Recall: " + str(round(pr_m.iloc[50]['Recall'], 2)))
    
    
    def meta_model_fitting(self):
        # Ignore warnings - need to find a better way of suppressing output later on
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        
        # Get button pushed on tab - tells us which ML fitting mode we're entering
        sender = self.sender()
        
        # CREATE MODELS + SAMPLING METHODS FOR LATER USE
        # create our models
        dummy = DummyClassifier(strategy='most_frequent')
        logit = LogisticRegression(max_iter=100)
        lda = LinearDiscriminantAnalysis()
        nb = ComplementNB()
        dtree = tree.DecisionTreeClassifier()
        knn = KNeighborsClassifier()
        svc = svm.SVC()
        bagged = BaggingClassifier()
        rand_forest = RandomForestClassifier()
        ext_trees = ExtraTreesClassifier()
        gboost = GradientBoostingClassifier()
        
        # create list of our models to loop through
        models = [dummy, logit, lda, nb, dtree, knn, svc, bagged, rand_forest, ext_trees, gboost]
        
        # create our sampling methods
        # oversamplers
        over = RandomOverSampler()
        smote = SMOTE()
        border_smote = BorderlineSMOTE()
        svm_smote = SVMSMOTE()
#        km_smote = KMeansSMOTE()
        adasyn = ADASYN()
        # undersamplers
        under = RandomUnderSampler()
        cnn = CondensedNearestNeighbour()
        tomek = TomekLinks()
        enn = EditedNearestNeighbours()
        n_cleaning = NeighbourhoodCleaningRule()
        onesided = OneSidedSelection()
        # combined samplers
        smoteenn = SMOTEENN()
        smotetomek = SMOTETomek()
        # none
        none = "none"
        
        # create list of sampling methods to loop through - removed ->  km_smote,
        samplings = [none, over, smote, border_smote, svm_smote,
                     adasyn, under, cnn, tomek, enn, n_cleaning,
                     onesided, smoteenn, smotetomek]
                     
                     
        # SET UP FOR META-LEARNING WITH METHYLSIGHT
        # We're on the meta ml fitting tab, meaning we need to alter our feature set
        print('altered feature set for meta learning.')   #db
        print('Onto generating our 2X training + 1X test (holdout) set for our meta learning model.')   #db
        
        # First get the selected basic model + sampling method combo
        b_model = self.base_model   # basic tuned model
        basic_samp = self.sampselect.currentIndex()   # sampling method
        
        if sender.text() == 'Fit Models and Balance Data - Return Best Fit':
            increment_by = 25/3
            QApplication.processEvents()
            self.auto_metafit_progress.setValue(round(increment_by))
            
        
        if samplings[basic_samp] == none:
            b_sampler = None
        else:
            b_sampler = samplings[basic_samp]
        
        print('initial length of x_features', len(self.feat_x))   #db
        print('initial length of y_features', len(self.feat_y))   #db
        
        # INCORPORATE METHYLSIGHT DATA WITH FEATURE SET DATA
        # First - identify holdout set from user data
        if 'SOURCE' in self.user_feat.columns:
            # User has included SOURCE column, indicating there is the presence of a holdout set
            self.m_holdout = self.user_feat[self.user_feat['SOURCE'] != 'PhosphoSitePlus']
            leftover = self.user_feat[self.user_feat['SOURCE'] == 'PhosphoSitePlus']
            
            # Create a new df of user-defined holdout values
            meta_test_x = pd.DataFrame(self.m_holdout[['SECONDARY_ML_SCORE']])
            holdout_x = pd.DataFrame(self.feat_x.loc[self.m_holdout.index])
            meta_test_y = pd.DataFrame(self.feat_y.loc[self.m_holdout.index])
            
            # Second - define training data for base and meta-learning models
            # first isolate features from leftover dataset
            temp_x = pd.DataFrame(self.feat_x.loc[leftover.index])
            temp_y = pd.DataFrame(self.feat_y.loc[leftover.index])
            
        else:
            # User has not included a SOURCE column, meaning they don't have a holdout set
            self.m_holdout = None
            
            # Create a holdout set from our full data
            temp_x, holdout_x, temp_y, meta_test_y = train_test_split(self.feat_x, self.feat_y, stratify=self.feat_y, test_size=0.3)
            meta_test_x = pd.DataFrame(self.user_feat.loc[holdout_x.index]['SECONDARY_ML_SCORE'])
            
        print('after removing holdout set - length of x_holdout', len(holdout_x))   #db
        print('after removing holdout set - length of y_holdout', len(meta_test_y))   #db
        print('after removing holdout set - length of temp_x', len(temp_x))   #db
        print('after removing holdout set - length of temp_y', len(temp_x))   #db
        
        
        temp_secondary_score = pd.DataFrame(self.user_feat.loc[temp_x.index]['SECONDARY_ML_SCORE'])
        temp_secondary_score = temp_secondary_score[~temp_secondary_score.index.duplicated(keep='first')]

        # Final split - split temp set into basic and meta ml training sets
        temp_x['SECONDARY_ML_SCORE'] = temp_secondary_score
        base_train_x, meta_train_x, base_train_y, meta_train_y = train_test_split(temp_x, temp_y, test_size = 0.5, stratify = temp_y)
        
        
        if sender.text() == 'Fit Models and Balance Data - Return Best Fit':
            QApplication.processEvents()
            self.auto_metafit_progress.setValue(round(increment_by*2))
        
        
        # Set up for fitting
        base_train_x = base_train_x.drop(columns=['SECONDARY_ML_SCORE'])  # Drop this - it isn't a feature for the base model
        
        meta_combo_train_x = pd.DataFrame(meta_train_x['SECONDARY_ML_SCORE'].tolist(), columns=['SECONDARY_ML_SCORE'], index = meta_train_x.index)   # Make a new df for meta learning
        meta_train_x = meta_train_x.drop(columns=['SECONDARY_ML_SCORE'])   # Can drop secondary score from meta now that it's safely in a new df
        
        
        # GENERATE SCORES FOR BASE ML MODEL USING SELECTED MODEL-SAMPLING COMBINATION
        # First fit model to base_train sets
        if b_sampler != None:
            base_train_x, base_train_y = b_sampler.fit_resample(base_train_x, base_train_y)
        base_model = b_model.fit(base_train_x, base_train_y)
        
        # Now predict scores using meta_train sets
        base_scores = base_model.predict_proba(meta_train_x)
        print(base_scores)
        base_scores = [item[1] for item in base_scores]
        meta_combo_train_x['PRIMARY_ML_SCORE'] = base_scores
        
        print('Successfully generated meta training scores using the base model.')   #db
        
        # Predict primary ml scores for holdout set
        base_holdout_scores = base_model.predict_proba(holdout_x)
        base_holdout_scores = [item[1] for item in base_holdout_scores]
        meta_test_x['PRIMARY_ML_SCORE'] = base_holdout_scores
        
        
        if sender.text() == 'Fit Models and Balance Data - Return Best Fit':
            QApplication.processEvents()
            self.auto_metafit_progress.setValue(round(increment_by*3))
        
        
        # NOW ONTO MODEL FITTING
        
        # AUTO-FITTING MODEL + SAMPLING METHOD
        if sender.text() == 'Fit Models and Balance Data - Return Best Fit':
        
            print('Onto the auto-model fitting.')   # db
            
            ## FIT VARIOUS MODELS AND SAMPLING METHODS - SELECT FOR THE BEST ONE :) ##
            
            f1 = []   # scoring metric
            completed = increment_by*3   # update progress bar
            increment_by = (100/4)/len(models)   # progress bar
            completed = completed + increment_by

            
            # First loop through models list and generate f1 scoring
            for m in models:
                print('Meta model fitting')   #db
                # We're on the metalearning model - manually determine our metric
                meta_model = m.fit(meta_combo_train_x, meta_train_y)   # fit to training data
                meta_model_pred = meta_model.predict(meta_test_x)    # generate scores using our test set
                metric = precision_score(meta_test_y, meta_model_pred)   # manually determine precision
                completed+=increment_by
                QApplication.processEvents()
                self.auto_metafit_progress.setValue(round(completed))
                f1.append(metric)

            
            # Return model with biggest f1 score
            model_max_f1 = max(f1)
            model_max = models[f1.index(max(f1))]

            
            # Secondly we go through sampling methods
            
            f1_score = 0
            f1_s = []
            increment_by = (100/4)/len(samplings)
            print('now onto sampling methods')   #db
            
            for s in samplings:
                if s == none:
                    metric = model_max_f1
                else:
                    print('Meta model fitting')   #db
                    # We're on the metalearning model - manually determine our metric
                    print('BALANCING METHOD', str(s))
                    print('META COMBO TRAIN')   #db
                    print(Counter(meta_train_y))   #db
                    mc_train_x, mc_train_y = s.fit_resample(meta_combo_train_x, meta_train_y)
                    meta_model_full = model_max.fit(mc_train_x, mc_train_y)   # fit model using pipeline
                    meta_model_full_pred = meta_model_full.predict(meta_test_x)   # predict using model on test set
                    metric = precision_score(meta_test_y, meta_model_full_pred)   # manually determine precision
                completed+=increment_by
                QApplication.processEvents()
                self.auto_metafit_progress.setValue(round(completed))
                f1_s.append(metric)

            
            # Identify sampling method with biggest f1 score
            samp_max_f1 = max(f1_s)
            samp_max = samplings[f1_s.index(max(f1_s))]
            
            # Make pipeline for sampling method + model to return
            if samp_max == none:
                m_sampler = None
            else:
                m_sampler = samp_max
            m_final_model = model_max
            
            # Return final f1 score
            f1_score = samp_max_f1
            
            # Return strings of model and sampling
            model_str = str(model_max)
            sampling_str = str(samp_max)
            
            model_index = f1.index(max(f1))
            sampling_index = f1_s.index(max(f1_s))
            
            print('best model', model_str)   #db
            print('best data balancing method', sampling_str)   #db
                        
            # Update drop down menus with selected model and sampling methods
            self.meta_modselect.setCurrentIndex(model_index)
            self.meta_sampselect.setCurrentIndex(sampling_index)
            print('Metalearning autofitting complete.')  #db
        
        
        # USER-SELECTED MODEL + SAMPLING METHOD
  
        elif sender.text() == 'Re-Run Fitting with New Selection':
            
            QApplication.processEvents()   # update progress bar
            self.meta_user_modfit_progress.setValue(15)
            
            i_u_mod = self.meta_modselect.currentIndex()   # model
            i_u_samp = self.meta_sampselect.currentIndex()   # sampling method
            
            # user-defined pipeline
            if samplings[i_u_samp] == none:
                m_sampler = None
            else:
                m_sampler = samplings[i_u_samp]
                meta_combo_train_x, meta_train_y = m_sampler.fit_resample(meta_combo_train_x, meta_train_y)
            
            m_final_model = models[i_u_mod]

            meta_model_full = m_final_model.fit(meta_combo_train_x, meta_train_y)   # fit model using pipeline
            meta_model_full_pred = meta_model_full.predict(meta_test_x)   # predict using model on test set
            f1_score = precision_score(meta_test_y, meta_model_full_pred)   # manually determine f1 score
                            
            # Return strings of model and sampling
            model_str = str(models[i_u_mod])
            sampling_str = str(samplings[i_u_samp])
            
            QApplication.processEvents()   # update progress bar
            self.meta_user_modfit_progress.setValue(30)
            
        # HYPERPARAMETER SEARCH FOR MODEL OF BEST FIT
        if m_sampler != None:
            # Apply sampling method to training data
            sampled_x, sampled_y = m_sampler.fit_resample(meta_combo_train_x, meta_train_y)
        else:
            sampled_x = meta_combo_train_x
            sampled_y = meta_train_y
        
        if self.meta_modselect.currentIndex() != 0:
            # Be sure we aren't using the dummy classifier (no hyperparameters for that)
            m_params = tuning(self.meta_modselect.currentIndex(), m_final_model, sampled_x, sampled_y)
        
            # Apply hyperparameters
            meta_model_full.set_params(**m_params)
            print(m_final_model.get_params())
        
        self.meta_model_final = meta_model_full
            
        # GENERATE PR AND ROC CURVE - META LEARNING MODEL
        
        threshold = 0.0
        t_rec = []
        t_prec = []
        t_spec = []
        t_sens = []
        thresh = []

        user_increment_by = 70/100
        user_completed = 30 + user_increment_by
        increment_by = (100/4)/100
        
        while threshold <= 1:
            # calculate probability on testing data
            y_proba = meta_model_full.predict_proba(meta_test_x)
            y_proba = y_proba[:, [1]]   #select the probability for the positive case only
            
            # select just the scores over our threshold
            cond = y_proba >= threshold
            # convert to 1 or 0 values
            y_pred = np.where((y_proba>=threshold), 1, y_proba)
            y_pred = np.where((y_pred<threshold), 0, y_pred)
            # calculate recall
            recall = recall_score(meta_test_y, y_pred)
            t_rec.append(recall)
            # calculate precision
            precision = precision_score(meta_test_y, y_pred)
            t_prec.append(precision)
            # calculate specificity
            tn, fp, fn, tp = confusion_matrix(meta_test_y, y_pred).ravel()
            specificity = tn/(tn+fp)
            t_spec.append(specificity)
            sensitivity = tp/(tp+fn)
            t_sens.append(sensitivity)

            # update threshold
            thresh.append(threshold)
            threshold += 0.01
            
            
            
            if sender.text() == 'Fit Models and Balance Data - Return Best Fit':
                QApplication.processEvents()
                completed+=increment_by
                self.auto_metafit_progress.setValue(round(completed))
                print('COMPLETED:', completed)
            elif sender.text() == 'Re-Run Fitting with New Selection':
                QApplication.processEvents()
                user_completed+=user_increment_by
                self.meta_user_modfit_progress.setValue(round(user_completed))
                print('USER COMPLETED', user_completed)
            
        fpr, tpr, _ = roc_curve(meta_test_y,  meta_model_full_pred)
        fpr = fpr.tolist()
        tpr = tpr.tolist()
        auc_m = roc_auc_score(meta_test_y,  meta_model_full_pred)
        
        pr_m = pd.DataFrame(list(zip(thresh, t_rec, t_spec, t_prec, t_sens)), columns=['Threshold'
                                                                                      , 'Recall'
                                                                                      , 'Specificity'
                                                                                      , 'Precision'
                                                                                      , 'Sensitivity'])
        roc_m = pd.DataFrame(list(zip(fpr, tpr)), columns = ['fpr', 'tpr'])
        
        
        # UPDATE OUTPUT + GENERATE PLOTS
            
        # Update f1 score to user
        self.meta_f1_score_label.setText("Current F1 Score: " + str(round(f1_score, 2)))
        
        # Call PR and ROC curve generation
        ms = ''.join(filter(str.isalnum, model_str))
        ss =''.join(filter(str.isalnum, sampling_str))
        save_annotation = 'runfiles/temp/Meta_ML_Model_' + ms + '_' + ss
        Plot(pr_m, (save_annotation+'_pr_curve'), roc_m, auc_m, (save_annotation+'_roc_curve'))
        
        # Update PR and ROC curves to user
        self.meta_pr_label.setPixmap(QPixmap(save_annotation+'_pr_curve.png'))
        self.meta_roc_label.setPixmap(QPixmap(save_annotation+'_roc_curve.png'))
        
        # Paths to images for file saving later
        self.m_pr = save_annotation+'_pr_curve.pdf'
        self.m_roc = save_annotation+'_roc_curve.pdf'
        
        # Update metric scores
        self.meta_sens_lab.setText("Sensitivity: " + str(round(pr_m.iloc[50]['Sensitivity'], 2)))
        self.meta_spec_lab.setText("Specificity: " + str(round(pr_m.iloc[50]['Specificity'], 2)))
        self.meta_prec_lab.setText("Precision: " + str(round(pr_m.iloc[50]['Precision'], 2)))
        self.meta_rec_lab.setText("Recall: " + str(round(pr_m.iloc[50]['Recall'], 2)))
        
        # Allow user to apply meta learning in next stage
        self.meta_model_exists = True
    
    def exp_proceed_upload(self):
        # USER SELECTED DATASET FROM LIST AND WANTS TO PROCEED (i.e. no metalearning unless it's methylsight for now)
        # Get which file they wanted to proceed with
        desired_exp = self.exp_dataset.currentIndex()
        if desired_exp == 0:
        # Lysines
            print("User selected to proceed with the surface exposed lysine dset")
            # upload it from its location
            lys_dset_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'datasets/dummy_stuff/dummyset_set8_features_025.csv')) # REPLACE LATER WITH FULL DSET!
            df = pd.read_csv(lys_dset_file)
        elif desired_exp == 1:
        # Arginine
            print("User selected to proceed with the surface exposed arginine dset")
            arg_dset_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'datasets/dummy_stuff/dummyset_set8_features_025.csv')) # REPLACE LATER WITH EXPOSED ARGININE DATASET!
            df = pd.read_csv(arg_dset_file)
            
        elif desired_exp == 2:
        # Tyrosine
            print("User selected to proceed with the surface exposed tyrosine dset")
            tyr_dset_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'datasets/dummy_stuff/dummyset_set8_features_025.csv')) # REPLACE LATER WITH EXPOSED TYROSINE DATASET!
            df = pd.read_csv(tyr_dset_file)
        elif desired_exp == 3:
        # Threonine
            print("User selected to proceed with the surface exposed threonine dset")
            thr_dset_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'datasets/dummy_stuff/dummyset_set8_features_025.csv')) # REPLACE LATER WITH EXPOSED THREONINE DATASET!
            df = pd.read_csv(thr_dset_file)
        self.user_exp_set = df
        # Run feature generation
        self.exp_feature_generate()
    
    def exp_upstate(self):
        # UPLOAD DATASET FROM USER
        upload_from = QFileDialog.getOpenFileName(self, "Upload Dataset to Run Through Machine Learning Model")
        upload_from_dir = upload_from[0]
        
        if ".xls" in upload_from_dir:
            # Excel file
            self.user_exp_set = pd.read_excel(upload_from_dir)
            self.user_experimental_label.setText('Uploaded excel file: ' + upload_from_dir)
            self.exp_feature_generate()
        elif ".csv" in upload_from_dir:
            # csv file
            self.user_exp_set = pd.read_csv(upload_from_dir, index_col=0)
            self.user_experimental_label.setText('Uploaded csv file: ' + upload_from_dir)
            self.exp_feature_generate()
        elif ".tsv" in upload_from_dir:
            self.user_exp_set = pd.read_csv(upload_from_dir, sep='\t', index_col=0)
            self.user_experimental_label.setText('Uploaded tsv file: ' + upload_from_dir)
            self.exp_feature_generate()
        else:
            # improper format >:(
            self.user_experimental_label.setText('Unable to upload file. Please attempt with a .xlsx, .csv, or .tsv file.')
            self.user_exp_set = None
    

    def exp_feature_generate(self):
        # DOWNSTATE FOR EXPERIMENTAL FEATURE GENERATE BUTTON - INITIATE FEATURE GENERATION SCRIPT FOR PASSED DATA FILE
        
        print("Running feature generation for experimental dataset")
        
        # Format user defined features for FeatureGen call
        sequences = self.user_exp_set['SITE_+/-7_AA']
        if self.user_exp_set['MOD_RSD'].dtype == 'int64':
            temp_annotation = self.user_exp_set['MOD_RSD'].astype(str)
        else:
            temp_annotation = self.user_exp_set['MOD_RSD'].str.split("-", expand=True)[0]
            temp_annotation = temp_annotation.str[1:]
        uid_annotation = self.user_exp_set['ACC_ID'] + "_" + temp_annotation
        
        self.user_exp_set = self.user_exp_set.set_index(uid_annotation)
        
        # Now get into generating our features!
        
        # Create df for results to go into
        v, h = FeatureGen(sequences[0], self.protdcal)
        features = pd.DataFrame(columns=h)
        features.loc[len(features)] = v
        
        increment_by = 100/len(sequences)
        completed = increment_by
        i = 1
        
        # Go through rest of sequences to generate feature set
        while i < len(sequences):
            ts = sequences[i]
            value, header = FeatureGen(ts, self.protdcal)
            features.loc[len(features)] = value
            completed+=increment_by
            i+=1
            QApplication.processEvents()
            self.exp_feat_progress.setValue(completed)
            self.exp_post_feat.setText("Dataset Feature Generation is " + str(round(completed, 1)) + "% complete.")

        # Make the index the same as our initial dataframe
        self.exp_set = features.set_index(uid_annotation)
        
        QApplication.processEvents()
        self.exp_feat_progress.setValue(100)
        
        self.exp_init_run.setEnabled(True)
        
        # Check for metalearning scores in experimental dataset
        if 'SECONDARY_ML_SCORE' not in self.user_exp_set.columns:
            self.exp_post_feat.setText("Feature Generation for the defined dataset has finished! \nThe chosen dataset has " + str(len(self.exp_set)) +  " features. You may proceed with Basic ML prediction.\nUnable to run Meta Learning prediction as the SECONDARY_ML_SCORE column in the uploaded dataset does not exist.\nTo proceed with Meta Learning, create the SECONDARY_ML_SCORE column with numeric scores from another computational method.")
            self.meta_model_exists = False
        elif (self.user_exp_set['SECONDARY_ML_SCORE'] == 0).all():
            self.exp_post_feat.setText("Feature Generation for the defined dataset has finished! \nThe chosen dataset has " + str(len(self.exp_set)) +  " features. You may proceed with Basic ML prediction.\nUnable to proceed with Meta Learning prediction as the SECONDARY_ML_SCORE column in the uploaded dataset has no values.\nTo proceed, populate the SECONDARY_ML_SCORE column with numeric scores from another computational method.")
            self.meta_model_exists = False
        else:
            self.exp_post_feat.setText("Feature Generation for the defined dataset has finished! \nThe chosen dataset has " + str(len(self.exp_set)) +  " features. You may proceed with Basic ML prediction.\nSECONDARY_ML_SCORE successfully detected. You may proceed with Meta Learning for the experimental dataset!")
            self.b2.setCheckable(True)   # user can now select meta learning model for experimental prediction run
        
        
    
    def experimental_ml_run (self):
        self.exp_init_run.setEnabled(False)   # false until we're done running this
        
        QApplication.processEvents()
        self.exp_ml_progress.setValue(10)
        
        # create our sampling methods
        # oversamplers
        over = RandomOverSampler()
        smote = SMOTE()
        border_smote = BorderlineSMOTE()
        svm_smote = SVMSMOTE()
        km_smote = KMeansSMOTE()
        adasyn = ADASYN()
        # undersamplers
        under = RandomUnderSampler()
        cnn = CondensedNearestNeighbour()
        tomek = TomekLinks()
        enn = EditedNearestNeighbours()
        n_cleaning = NeighbourhoodCleaningRule()
        onesided = OneSidedSelection()
        # combined samplers
        smoteenn = SMOTEENN()
        smotetomek = SMOTETomek()
        # none
        none = "none"
        
        QApplication.processEvents()
        self.exp_ml_progress.setValue(30)
        
        # create list of sampling methods to loop through - removed ->
        samplings = [none, over, smote, border_smote, svm_smote, km_smote,
                     adasyn, under, cnn, tomek, enn, n_cleaning,
                     onesided, smoteenn, smotetomek]
                     
        
        # Meta-Learning or Basic Model?
        if self.b1.isChecked():
            print('Basic Model')   #db
            # BASIC MODEL
            # Fit model to all training data
            # Get current model selection
            selected_mod = self.base_model # current model
            i_samp = self.sampselect.currentIndex()
            selected_samp = samplings[i_samp]   # current sampling method
            QApplication.processEvents()
            self.exp_ml_progress.setValue(40)
            
            
            # Fit to training data
            if selected_samp != "none":
                final_exp_samp = selected_samp
                exp_feat_x, exp_feat_y = final_exp_samp.fit_resample(self.feat_x, self.feat_y)
                model_full = selected_mod.fit(exp_feat_x, exp_feat_y)
            else:
                model_full = selected_mod.fit(self.feat_x, self.feat_y)
            
            
            QApplication.processEvents()
            self.exp_ml_progress.setValue(50)
            
            QApplication.processEvents()
            self.exp_ml_progress.setValue(60)
            
            # Predict on experimental data - add scores to file
            model_proba = model_full.predict_proba(self.exp_set)   # predict using model on test set
            QApplication.processEvents()
            self.exp_ml_progress.setValue(70)
            model_pred_y = model_proba[:, 1]  # isolate the positive scores only
            self.exp_set['PRIMARY_ML_SCORE'] = model_pred_y
            QApplication.processEvents()
            self.exp_ml_progress.setValue(80)
            
            # Sort by ML scores
            final_exp_scores = self.exp_set.sort_values(by=['PRIMARY_ML_SCORE'], ascending=False)
            QApplication.processEvents()
            self.exp_ml_progress.setValue(90)
            
            # Tell user how many predicted positives and negatives they have with a 0.5 cutoff
            positives = len(final_exp_scores[final_exp_scores['PRIMARY_ML_SCORE']>=0.5])
            negatives = len(final_exp_scores[final_exp_scores['PRIMARY_ML_SCORE']<0.5])
            QApplication.processEvents()
            self.exp_ml_progress.setValue(95)
            
        elif self.b2.isChecked():
            print('Meta Model')   #db
            # META MODEL
            # Format training data
            
            QApplication.processEvents()
            self.exp_ml_progress.setValue(40)
            
            # Define training data
            meta_train_full_x = self.feat_x
            meta_train_full_y = self.feat_y
            meta_train_full_x['SECONDARY_ML_SCORE'] = pd.DataFrame(self.user_feat[['SECONDARY_ML_SCORE']])
            
            print('length of full meta training set x-values', len(meta_train_full_x))
            print('length of full meta training set y-values', len(meta_train_full_y))
            print('length of full meta training set secondary scores', len(self.user_feat.loc[self.feat_x.index]['SECONDARY_ML_SCORE']))
            print('length of full meta training set secondary scores', len(self.user_feat[['SECONDARY_ML_SCORE']]))
            
            QApplication.processEvents()
            self.exp_ml_progress.setValue(50)
            
            # Split up training data
            meta_train_x, base_train_x, meta_train_y, base_train_y = train_test_split(meta_train_full_x, meta_train_full_y, test_size = 0.5, stratify=meta_train_full_y)
            
            # Get our basic model configuration
            selected_mod = self.base_model   # current model
            i_samp = self.sampselect.currentIndex()
            selected_samp = samplings[i_samp]   # current sampling method

            base_train_x = base_train_x.drop(columns=['SECONDARY_ML_SCORE'])
            
            if selected_samp == "none":
                basic_model = selected_mod.fit(base_train_x, base_train_y)
            else:
                base_train_x_s, base_train_y_s = selected_samp.fit_resample(base_train_x, base_train_y)
                basic_model = selected_mod.fit(base_train_x_s, base_train_y_s)
            
            
            QApplication.processEvents()
            self.exp_ml_progress.setValue(60)
            
            # Predict scores of 2/2 training data with basic model
            ms_meta_train_x = pd.DataFrame(meta_train_x[['SECONDARY_ML_SCORE']])
            meta_train_x = meta_train_x.drop(columns = ['SECONDARY_ML_SCORE'])
            basic_model_train_proba = basic_model.predict_proba(meta_train_x)   # predict using model on test set
            basic_meta_train_x = pd.DataFrame(basic_model_train_proba[:, 1] , columns=['PRIMARY_ML_SCORE']).set_index(meta_train_x.index) # isolate the positive scores only
            basic_meta_train_x['SECONDARY_ML_SCORE'] = ms_meta_train_x   #incorporate methylsight scores
            
            QApplication.processEvents()
            self.exp_ml_progress.setValue(70)
            
            # Predict scores of experimental data with basic model
            just_ms_exp = pd.DataFrame(self.user_exp_set[['SECONDARY_ML_SCORE']])
            basic_exp_proba = basic_model.predict_proba(self.exp_set)   # predict using model on experimental data
            basic_exp_x = pd.DataFrame(basic_exp_proba[:, 1] , columns=['PRIMARY_ML_SCORE']).set_index(self.exp_set.index) # isolate the positive scores only
            basic_exp_x['SECONDARY_ML_SCORE'] = just_ms_exp   #incorporate methylsight scores
            
            QApplication.processEvents()
            self.exp_ml_progress.setValue(75)
            
            # Get our meta model configuration
            m_selected_mod = self.meta_model_final   # current model
            m_samp = self.meta_sampselect.currentIndex()
            m_selected_samp = samplings[m_samp]   # current sampling method
            # Fit to training data
            if m_selected_samp == "none":
                meta_model = m_selected_mod.fit(basic_meta_train_x, meta_train_y)
            else:
                basic_meta_train_x_s, meta_train_y_s = m_selected_samp.fit_resample(basic_meta_train_x, meta_train_y)
                meta_model = m_selected_mod.fit(basic_meta_train_x_s, meta_train_y_s)
            
            QApplication.processEvents()
            self.exp_ml_progress.setValue(85)
            
            # Predict scores of experimental data with meta model
            meta_experimental_proba = meta_model.predict_proba(basic_exp_x)   # predict using model on experimental set
            meta_experimental_scores = meta_experimental_proba[:, 1]
            
            QApplication.processEvents()
            self.exp_ml_progress.setValue(90)
            
            # Add scores to file
            temp_exp = self.user_exp_set
            temp_exp['PRIMARY_ML_SCORE'] = basic_exp_x['PRIMARY_ML_SCORE']
            temp_exp['META_LEARNING_SCORE'] = meta_experimental_scores
            
            # Sort from biggest to lowest
            final_exp_scores = temp_exp.sort_values(by=['META_LEARNING_SCORE'], ascending=False)
            
            QApplication.processEvents()
            self.exp_ml_progress.setValue(93)

            # Tell user how many predicted positives and negatives they have with a 0.5 cutoff
            positives = len(final_exp_scores[final_exp_scores['META_LEARNING_SCORE']>=0.5])
            negatives = len(final_exp_scores[final_exp_scores['META_LEARNING_SCORE']<0.5])
            QApplication.processEvents()
            self.exp_ml_progress.setValue(95)
        
        self.exp_post_run.setText('Successfully completed prediction.\nWith a threshold of 0.5, there are ' + str(positives) + ' predicted positives and ' + str(negatives) + ' predicted negatives.')
        
        # Progress bar for ml run
        QApplication.processEvents()
        self.exp_ml_progress.setValue(100)
        
        self.experimental_scores = final_exp_scores
        
        self.exp_init_run.setEnabled(True)   # done run so we set back to true
        
        self.exp_save.setEnabled(True)   # done run so we can now save the score file :)
        
        
    def experimental_ml_run_save(self):
        exp_scores_dest = QFileDialog.getExistingDirectory(self, 'Save Score File To')
        filesavetime = time.strftime("%Y-%m-%d_%H.%M.%S")
        self.experimental_scores.to_csv(exp_scores_dest+'/ProteoML_Genius_predictions_' + filesavetime + '.csv')
    
    
def FeatureGen (sequence, protdcal):
    # METHOD FOR GENERATING FEATURES USING ONE SEQUENCE AT A TIME
    # FIRST: GENERATE PROTDCAL VALUES
    slist = list(sequence)   # first split sequence up into list
    # Go through sequence to get protdcal value
    pd = []
    for i in slist:
        pd.append(protdcal.loc[i].tolist())
    values = list(map(lambda *x: sum(x), *pd))   # add up values
    headers =  protdcal.columns.tolist()   # include headers
    
    
    # SECOND: GENERATE ONE-HOT ENCODING
    aa = ['K', 'R', 'H', 'A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y',
          'N', 'C', 'Q', 'S', 'T', 'D', 'E', 'G', 'P']   # possible amino acids
    # Make headers and one-hot encoding for each letter
    for i in aa:
        j = 0
        while j < len(sequence):
            headers.append('ONE-HOT_' + str(j) + '-' + i)  # make header
            if sequence[j] == i:
                values.append(1)
            else:
                values.append(0)
            j+=1
    
    
    # THIRD: GENERATE MACCS KEYS
    # Generate maccs keys
    mol = (rdmolfiles.MolFromFASTA(sequence))
    fp = (MACCSkeys.GenMACCSKeys(mol))
    maccs = fp.ToBitString()
    binary = list(maccs)   # split up into list
    values.extend(binary)   # add list onto resulting values
    # Generate headers for maccs keys
    mt = list(itertools.chain(range(len(binary))))
    mt = [str(s) + '_maccs' for s in mt]
    headers.extend(mt)   # append header values
    
    return values, headers

def k_fold(X, y, model, sampling, sampling_type):
    # K Fold cross validation for model fitting - returns f1 score as metric

    if sampling != 'none':
        steps = [(sampling_type, sampling), ('model', model)]
        pipeline = Pipeline(steps=steps)
        cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=1)
        scores = cross_val_score(pipeline, X, y, scoring='f1', cv=cv, n_jobs=1)
    else :
        cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=1)
        scores = cross_val_score(model, X, y, scoring='f1', cv=cv, n_jobs=1)
    score = mean(scores)
    
    return score
    
def tuning(model_not, model, sampled_x, sampled_y):
    # Tune our chosen model using a predefined space of hyperparameters dependent on the model type
    # Tuning will consider our sampling/balancing method as applied to our dataset
    
    print(model)
    
    # Big if/else for our model
    # Nothing for dummy classifier
    
    if model_not == 1:
        # Logistic Regression
        parameters = {'penalty':('l2', None), 'solver':('lbfgs', 'sag', 'liblinear'), 'max_iter':[10, 50, 100, 1000]}
        search = 'random'
        
    elif model_not == 2:
        # Linear Discriminant Analysis
        parameters = {'solver':('svd', 'lsqr', 'eigen')}
        search = 'grid'
        
    elif model_not == 3:
        # Complement Naive Bayes
        parameters = {'alpha':[1e-5, 1e-3, 0.1, 0.5, 1], 'norm':(True, False)}
        search = 'grid'
        
    elif model_not == 4:
        # Decision Tree Classifier
        parameters = {'criterion':('gini', 'entropy', 'log_loss'), 'splitter':('best', 'random'), 'min_samples_split':[2, 5, 10, 20, 40], 'min_samples_leaf':[1, 2, 4, 8, 16, 20]}
        search = 'random'
        
    elif model_not == 5:
        # K-Nearest Neighbours Classifier (k-NN)
        parameters = {'weights':('uniform', 'distance', None), 'leaf_size':list(range(1,50)), 'n_neighbors':list(range(1,30)), 'p':[1, 2]}
        search = 'random'
        
    elif model_not == 6:
        # Support Vector Classifier (SVC)
        parameters = {'C':[0.01, 0.1, 1, 10], 'gamma':('scale', 'auto'), 'kernel':('linear', 'rbf', 'poly')}
        search = 'random'
        
    elif model_not == 7:
        # Bagging Classifier
        parameters = {'n_estimators':[1, 2, 4, 8, 16]}
        search = 'grid'
        
    elif model_not == 8:
        # Random Forest Classifier
        parameters = {'n_estimators':[100, 200, 500, 1000, 2000], 'max_depth':[10, 20, 50, 100], 'min_samples_split':[2, 5, 10], 'min_samples_leaf':[1, 2, 4], 'max_features':('sqrt', 'log2'), 'bootstrap':(True, False)}
        search = 'random'

    elif model_not == 9:
        # Extra Trees Classifier
        parameters = {'n_estimators':[100, 200, 500, 1000, 2000], 'min_samples_leaf':[5, 10, 20], 'max_features':[2, 3, 4]}
        search = 'random'
        
    elif model_not == 10:
        # Gradient Boosting Classifier
        parameters = {'max_depth':[3, 5, 7, 9, 10], 'n_estimators':[1, 2, 5, 10, 20, 50, 100, 200, 500], 'learning_rate':loguniform.rvs(0.01, 1, size=10).tolist()}
        search = 'random'
        
    if search == 'random':
        # Set up randomsearchCV
        clf = RandomizedSearchCV(model, parameters)
    else:
        # Set up gridsearchCV
        clf = GridSearchCV(model, parameters)
    # run the search
    result = clf.fit(sampled_x, sampled_y)
    print(result.best_params_)
    # return dict of best hyperparameters for our model
    return result.best_params_

    

def Plot (pr_metrics, pr_saveas, roc_metrics, roc_auc, roc_saveas):
    pr_metrics = pr_metrics.drop(columns='Sensitivity')
    # Plot scores from metrics method as pr and roc curves
    sns.set_context('paper')
    ax = sns.lineplot(x='Threshold', y='value', hue='variable', data=pd.melt(pr_metrics, 'Threshold'), ci=None)
    ax.set(ylabel='Performance')
    ax.grid()
    ax.legend(title='', bbox_to_anchor=(.5, 1), loc='lower center', ncol=3)
    plt.ylim(0,1)
    plt.xlim(0,1)
    pr_save_pdf = os.path.abspath(os.path.join(os.path.dirname(__file__), (pr_saveas+'.pdf')))
    pr_save_png = os.path.abspath(os.path.join(os.path.dirname(__file__), (pr_saveas+'.png')))
    plt.savefig(pr_save_png, dpi=72, bbox_inches='tight')
    plt.savefig(pr_save_pdf, dpi=300, bbox_inches='tight')
    plt.clf()
    plt.plot(roc_metrics['fpr'],roc_metrics['tpr'],label="Model, auc="+str(roc_auc))
    plt.legend(title='', bbox_to_anchor=(.5, 1), loc='lower center', ncol=1)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    roc_save_pdf = os.path.abspath(os.path.join(os.path.dirname(__file__), (roc_saveas+'.pdf')))
    roc_save_png = os.path.abspath(os.path.join(os.path.dirname(__file__), (roc_saveas+'.png')))
    plt.savefig(roc_save_png, dpi=72, bbox_inches='tight')
    plt.savefig(roc_save_pdf, dpi=300, bbox_inches='tight')
    plt.clf()

if __name__ == '__main__':
    app = QApplication([])
    window = Window()
    window.show()
    sys.exit(app.exec())
