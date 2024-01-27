# SoloKubayashi-Maru-Prostate-Cancer-Patients-

My name is Leen Madani and I am a student enrolled in the Master of Biotechnology program (Digital Health Technologies stream) at the University of Toronto.

My repositories contain projects that I worked on throughout my master's year. If the project is done with teams, then this will be clearly indicated.

The projects are primarily relevant to data science in health, where I attempt to answer a research question through a thorough analysis of a dataset. The process involves data cleaning and the usage of a lot of statistical techniques (regression models, survival analysis, machine learning, classification, etc..).

If you have any questions, please email me at leen.madani@mail.utoronto.ca

********************
Solo Project & Presentation 

The objective of your analysis is to describe any potential effect of different treatments to health utility, in prostate cancer patients. In this longitudinal study data from prostate cancer patients (found here Download here) have been collected, at three time points. a) At T1, baseline, prior to treatment, b) at T2, 3 months after treatment and c) at T3, 1 year after treatment. The data include utility scores based on different kinds of instruments (EQ5D, HUI2, HUI3 and PORPUS-U), co-morbidities, disease progression related information, and type of treatment. You will focus on EQ5D as a measure of health utility. Some of the variables included in this data set and are of interest in this analysis are:

 

Variable

Description

AGE

Age in years

COHORT

Cohort (A=newly diagnosed; B=metastatic; C=follow-up)

DXSTAGE

Diagnosis stage

T1CHEMO

Received chemotherapy (0=no; 1=yes; 2=not sure)

T1DIABET

Have diabetes (0=no; 1=yes (past, not now); 2=Yes (now); 3=not sure)

T1EQTOT

EQ-5D utility (range=0-1; 1 best)

T1KIDNEY

Have kidney problems (0=no; 1=yes (past, not now); 2=Yes (now); 3=not sure)

T1RP

Had radical prostatectomy (0=no; 1=yes; 2=not sure)

T1RT

Had radiation therapy (0=no; 1=yes; 2=not sure)

T1SPBONE

Has tumour spread to bones 0=no; 1=yes; 2=not sure

T1VIAGRA

Taken Viagra 0=no; 1=yes; 2=not sure

ptid

PATIENT ID

t1ADT

Taken any type of hormone (ADT) 0=no;1=yes

t1arthritis

Have disabling arthritis (0=no; 1=yes (past, not now); 2=Yes (now); 3=not sure)

t1heart

Have heart problems (0=no; 1=yes (past, not now); 2=Yes (now); 3=not sure)

t1sporg

Has tumour spread to other organs 0=no; 1=yes; 2=not sure

t2ADT

Taken any type of hormone (ADT) 0=no;1=yes

t2arthritis

Have disabling arthritis (0=no; 1=yes (past, not now); 2=Yes (now); 3=not sure)

t2rp

Had radical prostatectomy (0=no; 1=yes; 2=not sure)

t2rt

Had radiation therapy (0=no; 1=yes; 2=not sure)

t3ADT

Taken any type of hormone (ADT) 0=no;1=yes

t3arthritis

Have disabling arthritis (0=no; 1=yes (past, not now); 2=Yes (now); 3=not sure)

t3rp

Had radical prostatectomy (0=no; 1=yes; 2=not sure)

t3rt

Had radiation therapy (0=no; 1=yes; 2=not sure)

 

The main treatments are radical prostatectomy (RP) and radiation therapy (RT). For each one of those, hormone therapy (ADT) could have been followed as well as a secondary treatment.

With the exception of age, ptid, cohort and dxstage, all other variables are measured in 3 time points. If a variable is not included in the above table, it is assumed that its name has the same format as the ones that are included. E.g. in addition to T1EQTOT, in the data there are also T2EQTOT and T3EQTOT. You need to focus only on cohort A, newly diagnosed patients, therefore all other patients need to be ignored in the analysis.

Choose the appropriate methods for approaching the problem and for analyzing the data, in order to answer the research question in a complete and comprehensive way. You will present your results in a short (concise and informative) 5 min presentation of 5 slides maximum.
