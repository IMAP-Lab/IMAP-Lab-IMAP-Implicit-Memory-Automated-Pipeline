Written by R.J. Addante, 3/26/2024

Note on event codes for data used in Muller, Sirianni, & Addante, 2021; Addante et al., 2023; Addante et al., 2024: 
-First: be aware that extensive details are provided in the IMAP script files
-Second: Experimental structure of the memory retrieval paradigm had event codes logged in the EEG file (and behavioral logfiles, too) for the stimuli onset, as well as for the responses provided by participants. Details for this are expanded upon below.

Stimulus onset events: if the item (a word) was an old item from the study phase, and was from the animacy task (see Addante et al., 2023, Psychophysiology; Muller, Sirianni, & Addante, 2021, European Journal of Neuroscience), then its stimulus onset was logged with a code of "61". If the item was old (from the study/encoding phase) and was instead from the other task ('Manmade' decision task), then it was logged with a 62. If the item was a new word (not presented earlier in the encoding phase) it was logged with a 63. Thus:
---> 61 = Old item, animacy task
---> 62 = Old item, manmade task
---> 63 = New item, no study task

Memory responses were then made for each item: 
A) an item recognition confidence judgment on a scale of 1 to 5, with 1 being 'sure new', and 5 meaning 'sure old'
_________________________________________________________________
New                                                 Old

1        2                3                4         5

Sure    Unsure                            Unsure    Sure
_________________________________________________________________

This response above was then followed immediately by the response/promt below:

B) a source memory confidence task that was also 5-points of confidence. 

For the source memory confidence task, the response options/values were as follows: 
_________________________________________________________________
Animacy Task        Source Unknown        Manmade Task

1        2                3                4         5

Sure    Unsure      Source Unknown        Unsure    Sure
_________________________________________________________________

Thus: event code conditions that start with 61 (Animacy condition) will have correct source judgment response code of '1' and '2', and incorrect source judgments of '4' and '5', whereas the opposite will be true for the event code condition of 62 (manmade task) in which correct source responses will be '4' and '5' and incorrect source responses will be noted as '1' and '2'
(note the contrasting scale for encoding task conditions described above, e.g. in *carefully* writing binoperation files in ERPLab)

Then: 
As noted in Muller et al (2021, EJN): every ten retrieval trials there was also a metacognitive judgment asked for particpants, e.g. what percentile do you think you are performing in on the  memory test?  
Subjects responded on a 1 to 5 confidence scale indicating which percentile they thought they were performing in:
_________________________________________________________________
< 60's    60's     70's        80's      90's

1          2        3           4         5
_________________________________________________________________
For this trial's stimulus onset, the event code was an 81, and the response given by the participant (1 through 5 rating) was logged as the ensuing event code. 

So, overall, sequences of event codes could look as follows:
-Stimulus onset (61/62/63)
-Item recognition confidence: 1/2/3/4/5
-Source memory confidence: 1/2/3/4/5 (note the contrasting scale for encoding task conditions described above, e.g. in *carefully* writing binoperation files in ERPLab)
-Metacognition Prompt: 81
-Metacognitive Resposne: 1/2/3/4/5
-----> Sample code series of 61, 1, 5, 81, 4 [6115814]

Note on EEG Data, 1: as noted in the IMAP scripts' comments (though also see publication Methods of Addante et al., 2023, Psychophysiology; Muller, Sirianni, & Addante, 2021, European Journal of Neuroscience for more specific technical details): the EEG data has been preprocessed as being re-referenced offline to the average of the mastoids, High pass filtered at .1 Hz, downsampled to 256 Hz, applied channel operations, and with events listed in ERPab,. and epoched. This EEG data has also had independent components analysis (ICA) performed on it, followed first by artifact correction, and then by manual artifact rejection for any remaining artifacts identified by researchers in a manner unbiased to experimental/memory conditions (e.g. blind). This EEG data is ready to then be created into ERP files (i.e. in ERPLab, or comparable program), apply binoperations (in ERPLab, for instance), etc. for post-processing. 
Note on EEG Data, 2: as noted in aforementioned publications, several participants were excluded from analysis due to various technical reasons (Addante et al., 2023), so they remain excluded from the uploaded data.
