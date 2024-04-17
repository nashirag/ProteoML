# ProteoML
A desktop application to uncover post-translational modifications specific to a particular enzyme via automated machine learning.

ProteoML is available [here for mac](https://drive.google.com/file/d/14z-RbkTZDaD6EShj9_ObrBQwLnWY5Lpu/view?usp=share_link) and [here for PC!](https://drive.google.com/file/d/1R_sth_UdITcKrPxIyMX80K4e0fx0iPGv/view?usp=sharing)

You can find the accompanying application, [ProteoML Visualizer here.](https://drive.google.com/file/d/1KF31HU1NzGKHyO8b_0nECHdwLd5cHAoM/view?usp=share_link)

## Installation Guide - Mac OSX (tested on Sonoma 14.2 and up)
[Download the compressed ProteoML application file here.](https://drive.google.com/file/d/14z-RbkTZDaD6EShj9_ObrBQwLnWY5Lpu/view?usp=share_link)

After unzipping the file, move it to your applications folder, then open it. It is likely that OSX security permissions will not allow you to open the file without some extra steps, as found [here](https://support.apple.com/en-ca/102445).

_NOTE: ProteoML will not damage your computer, and contains no malware. This alert is due to the refusal of all applications created by those without a developer account (a fee-based program)._

To subvert this warning, go to the Privacy & Security section of System Settings. For one hour after your attempt to open ProteoML, you will see a new option under the Security option stating that _'"ProteoML" was blocked from use because it is not from an identified developer.'_ Click "Open Anyway" to initiate a secondary popup stating that _'"ProteoML" can't be opened because Apple cannot check it for malicious software.'_ Click "Open" and _voila_!

If you instead get a popup stating that ProteoML is damaged, delete the application and attempt the install once again. [OSX will deem any application as "damaged" if the security permissions fail to be appropriately handled](https://support.apple.com/en-ca/102445). 


## Installation Guide - Windows 32-bit (tested on Windows 10 64-bit)
[Download the ProteoML installer here.](https://drive.google.com/file/d/1R_sth_UdITcKrPxIyMX80K4e0fx0iPGv/view?usp=sharing)

You will get a popup with the header _"Windows protected your PC"_ - click the underlined _"More Info"_ to allow a secondary _"Run anyway"_ button. Click this button to initiate the installer window - select _"Yes"_ once it asks you _"Do you want to allow this app from an unknown publisher to make changes to your device?"_ This will initiate the installer wizard. 

Follow the steps to install ProteoML locally, and _voila!_

##

## 1 – Download Peptide List
This tab allows the user to select and download files specific to a PTM of interest. 

Simply click to select the dataset(s) you want. Once selected, click Download File(s) to download the excel files locally to your machine. 

##

## 2 – Upload Peptide Results
Here users upload their experimentally tested dataset and generate features. 

After experimentally determining peptide sites that obtain PTMs from your enzyme of interest, upload an excel, .csv, or .tsv file by clicking the “Upload Result Dataset” button. This will automatically load the table and headers to ProteoML. You must then define which column holds the peptide sequence, protein, ID, experimental binary classification, etc. Click the “Initiate Feature Generation” button to start the feature generation process. 

**NOTE:** You must specify at least the Protein ID, AA Sequence, and Experimental fields with the corresponding columns to run feature generation. 
- _Protein ID:_ Uniprot ID, gene name, given name, etc. to track and differentiate sequences (Uniprot ID is ideal for later ProteoML Visualizer analysis)
- _AA Sequence:_ Amino acid sequence (i.e. GGARWLKPRTSGA, note central AA is modified AA)
- _Experimental:_ Binary score indicating enzymatic activity for the peptide tested (1 = activity, 0 = no activity)

**NOTE:** This process may take a while to run, depending on the size of your experimental dataset. An option is to set this up and let it run in the background, or overnight (depending on the size of your feature set!). A dataset of approximately 4,500 sites takes about 1 hour to run.

**NOTE:** To proceed with Meta Learning (tab 4), you must specify a Secondary Score column, containing the scores for each peptide/site within the protein determined by another ML scorer (i.e. MethylSight for kme, or MuSite Deep for 12 types of PTMs). The Validation score is for defining a holdout set specific to the secondary ML predictor. This will likely be unavailable for most secondary ML predictors, but it will ensure that when testing your model you aren’t potentially introducing bias. 

**NOTE:** If you mis-enter information here (i.e. mix up column headings), or your data is not of the correct format, the program may give you an error message, or exit prematurely. Please ensure that your data somewhat resembles that in Figure 1 (headings can vary, and only Protein ID, AA Sequence, and Experimental columns absolutely required). 

![Screenshot of a table.](/assets/scrngrb_file.png)

**Figure 1.** Example dataset for upload into ProteoML. UNI_ID column contains Uniprot ID, SITE_LOC contains the site location of the central lysine within the full length protein, SITE_+/-7_AA depicts the peptide sequence representative of the site, and EXPERIMENTALLY_ACTIVE holds the binary values of the experimental results for the enzyme tested. Secondary ML predictor scores (i.e. Musite Deep, MethylSight) are listed within the SECONDARY_ML_SCORE column, and DSET contains the distinction for the validation set to test the ensemble ML model with. 

##

## 3 – Fit Models
This tab is where the primary, or base ML model is fit to the experimental dataset provided. 

Clicking the “Fit Models and Balance Data – Return Best Fit” button initiates the auto ML functionality of ProteoML. This process will run through all available ML models and data balancing methods to identify the best-scoring combination for your dataset, based off of F-score for imbalanced datasets (i.e. your experimental data had a lot more negative sites of enzymatic activity than positive sites), or the ROC-AUC (Receiver Operating Characteristic – Area Under the Curve) for balanced datasets (i.e. experimental data showed relatively equal amounts of positive and negative data). After the best model and sampling combination is determined, the score will be output, along with a Precision-Recall Curve, and Threshold v Performance Metrics so you may determine the fit of your model. The “Metrics at Specified Threshold” is automatically updated to the optimal decision threshold for your model, however you can increase or decrease it if you’d prefer. 

If you’d like to select your own model, or alter the data balancing method or model type after an auto ML run, you may do so using the dropdown menus. The “Machine Learning Model” dropdown lists all possible models, and the “Data Balancing Method” lists all potential data balancing techniques you may apply within ProteoML (limited for the problem at hand). Clicking “Re-Run Fitting with New Selection” will re-run the metric and graphing generation for your new selection, also supplying a new cutoff threshold if necessary. 

The “Highly Associated GO Terms” table lists the molecular function GO terms for your training data, if you supplied a Uniprot ID as the Protein ID when you uploaded peptide results in tab 2. This tells you the most associated molecular functions with your dataset. 

**NOTE:** The model most recently selected for via the auto ML function, or the manual ML and data balancing function is the model used for the Ensemble ML model, and for the Experimental Dataset Scoring. If you want to change the model used in these processes, you will need to either click the “Fit Models and Balance Data – Return Best Fit” button or the “Re-Run Fitting with New Selection” buttons to ensure the new model is initialized before carrying on. This also applies to the selected threshold. If you wish to change the threshold, you must click “Update Metrics” after entering the new value to ensure this new threshold is applied. 

**NOTE:** Model fitting can take a long time, and some steps can seemingly “freeze” the program. Don’t worry if a run takes 30+ minutes, as the program is running complex ML fitting, scoring, and graphing processes in the background. 

##

## 4 – Meta Learning
This tab allows users to generate Ensemble or Meta Learning ML models. 

Much like the Fit Models tab, this tab offers similar functionalities (auto ML process and user selection process). However, instead of building an ML model from the training data, this process combines the scores generated from the selected base ML model optimized in tab 3, as well as the secondary ML scores defined within the uploaded training dataset within a new ML model. You’ll notice a new option is listed within the Machine Learning Model dropdown selection, “Soft Voting.” This process takes either the balanced or weighted mean of the ML scores of the base and secondary models. You may instead apply a new ML model rather than Soft Voting, however it is very important that complex ML models be avoided at this stage. Meta ML is a great way to reduce overfitting by incorporating more scorers, however if the model used to combine the scores is too smart, it will be overfit to your data. This means that your ensemble model will be very good at predicting the scores of your training data, but bad at predicting the scores of data it’s never seen before (testing/validation data, or experimental datasets later). 

Overfitting is avoided within the auto ML step by limiting the potential models to only the simplest: soft voting, logistic regression, and linear discriminant analysis. Within the Precision-Recall curves shown, there are two lines: one for model training, one for model testing. The model training line depicts the Precision-Recall performance of the ML model on the training dataset (i.e. data the model has already seen before), whereas the model testing line shows the performance for the testing/validation dataset (i.e. data the model has not yet seen). The distance between these two lines tells us the efficacy of our model. Generally, we expect to see a slight decrease in performance between from the training dataset to the testing dataset. If the testing dataset greatly decreases in performance from the training dataset, then our model is overfit. On the other hand, if our training dataset greatly resembles the performance of our testing dataset, then the model is underfit. The purpose of including both lines is so you can decide how well-fit the model is prior to scoring the experimental set. 

If you think there is too large of a gap between the performances of the model, try decreasing the intricacy of the ensemble model, and not applying sampling. If that does not improve the fit of the model, try going back to the Fit Models tab (tab 3), and reduce the complexity of the base model, or remove the data balancing altogether (select “None” from the dropdown menu). If the issue persists, it may be due to your secondary ML scorer. You may just use the base ML model (assuming that is fit correctly to your data) or try to get the ensemble model to the best state you can, then proceed.   

**NOTE:** This process is generally shorter, as just two features need to be fit within the ML models. 

**NOTE:** You will be unable to run this stage if you did not specify the secondary ML scoring column within the “Secondary Score” dropdown in tab 2 (Upload Peptide Results).

##

## 5 – Experimental Predictions
Score an “experimental” dataset for enzymatic activity using the constructed ML model. 

The first option is to employ the use of a pre-loaded dataset (contains MuSite Deep scores as secondary ML scores). You may download the selected dataset to alter the set itself, or include a different method of secondary ML scoring, or you may simply proceed with your selection (Proceed with Selected Dataset). If you choose to download this dataset and alter it in some way, you’ll need to reupload it under Option 2, by clicking the “Upload or Reupload a Dataset” button. This applies as well if you’re uploading a unique dataset to score. As with the initial feature generation step, you must specify the columns that contain the relevant data (Protein ID, Site, Peptide Sequence, and Secondary Score). 

Regardless of if you’ve selected to use a pre-loaded dataset, reupload, or upload your own, you must click the “Initiate Feature Generation for Experimental Set” button. Clicking this button will kick off another round of feature generation, which may take even longer than the first round for the training dataset (tab 2). Once this process is complete, you may select the “Basic Machine Learning Model” checkbox, or the “Meta Machine Learning Model” checkbox if you’ve initialized these models. Next click “Initiate Machine Learning Run” to score your file. Once complete (this process runs fast), click “Save Score File Locally” to save the .csv file in your specified directory. 

**NOTE:** If you’ve supplied ProteoML with Secondary Scores in the training dataset, as well as the experimental dataset (indicate the column via the Secondary Score dropdown menu), then you may select the “Meta Machine Learning Model” box. When your scores are output by ProteoML, both the primary/base ML model and secondary ML model scores will be output together with the ensemble ML model score. There is no need to re-run the scoring with the Basic Machine Learning Model to retrieve those scores.

**NOTE:** Regardless of the cutoff score supplied, the primary and/or secondary and ensemble ML scores will be output. The cutoff score defines the Final ML Classification column output in the save file (greater than cutoff = 1, less than = 0). 

**NOTE:** A pre-uploaded dataset exists of which features have been pre-generated— Surface Exposed Sites – Lysine (MethylSight kme Scoring) under the Select a Dataset dropdown menu. This dataset contains the surface exposed lysine proteome, however it is limited to peptides 15 AA in length. Other datasets with variable lengths will be populated soon. 

##

## Saving your Session
In the current state of ProteoML, you may only save your session by attempting to close the window or quit the application. File saving, and feature set re-uploading capability will be expanded upon shortly in the next iteration of the application. For now, to save the session, go to close the window and a popup will appear, asking if you’d like to save your progress. Click “Yes” and another popup will appear allowing you to specify where the file is saved to. Select the appropriate file on your local machine, and a file named “ProteoML_Genius_session...” will be saved there, with the time and date. 

The save file contains the following:

- _Feature importance graphs_ for your base ML model—these are features deemed to be important by your base ML model, as determined by permuted feature importance.

- _A logfile_ containing the logs of your session. This is for debugging purposes.

- _Performance metric curves_ for your base ML model and the ensemble model (if initialized). The precision-recall (PR) curve, receiver operating characteristic (ROC) curve, and the performance metric v. threshold curve will be listed for the models generated.

- _Feature sets_ for the training dataset and experimental sites/peptides for potential future analysis.
Result dataset containing the last scored experimental file by the ML model.


**NOTE:** Your ROC curves may look bad (i.e. close to the x=y line). That is to be expected, as the ROC curve does not consider imbalanced data. The performance of the ROC curve decreases steadily with less positive data. 


**NOTE:** A log file should be saved locally, save these if you encounter issues for debugging.

##

# ProteoML Visualizer

This application was created to visualize ProteoML predictions, or any other collection of proteins within a network of protein-protein interactions. 

Currently, the only available network to map to is the HuRI interactome, however there are plans to include the STRING-db interactome shortly, as well as to allow users to upload their own. In the future, the dropdown will be editable, but for now we maintain use of the HuRI interactome. 

Users select their list for download by selecting the “Click to upload a file containing Uniprot IDs” button, which will initiate the file upload dialogue. Once selected, the columns of the uploaded file (.xlsx, .csv, or .tsv are the accepted file types) will be added to the dropdown list to the right of “Column Containing Uniprot IDs”. As the title suggests, the user selects the column from their tabular data which contains the Uniprot ID alone. 

After selection, the user hits “Click to initiate SAFE run” to run the SAFE analysis for the uploaded Uniprot IDs. A network diagram will appear, and the “GO Terms Mapped to Network” table will become populated with the corresponding information. The user may save these both individually.

**NOTE:** Currently there is no readout of SAFE activity. This feature is under development and will be supported shortly. For now, it may seem as if the program is “frozen”, however that only signifies that SAFE is running. Leave the program to complete the run. This may take several hours, or overnight due to the complexities of the clustering networks applied..

**NOTE: We strongly suggest users upload positive sites of predicted enzymatic activity output by ProteoML as the data for the Visualizer.** This will help inform the accuracy of the ProteoML models. A secondary experiment may be to run SAFE analysis on the negatively predicted sites for comparison.

