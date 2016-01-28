# Introduction
**DatasetManager** is a singleton developed to retrieve dataset information for 
samples used by the IFCA and UO groups in CMS. This information is kept in a 
Google Spreadsheet that is regularly updated by our operators.

The information available for each sample currently is:
+ **cross section**
+ **number of events** in the sample
+ **list of files** composing the sample with full path. DatasetManager *guesses*
the path to the files based on the domain name of the machine where it is been run. 

DatasetManager also works for data samples where only the last bullet above is
relevant. By searching the path where files should be stored it is able to find
the files that compose our data samples.

# Instructions
## Getting the DatasetManager
As a singleton its **instance** may be retrieved by calling:

```cpp
DatasetManager* dm = DatasetManager::GetInstance();
```

## Selecting a tab in the spreadsheet
To **select** one of the many **tab**s in the spreadsheet use:

```cpp
dm->SetTab("Tab_name");
```
The first time this is done for a given tab, the information on the datasets is
**automatically downloaded** to a local cache (and stored in folder ```Datasets```).
Next time the local cache will be used to speed access. If the local cache is
older than 24 hours the information will be redownloaded. Downloading the information
from the web may always be forced by calling ```DatasetManager::RedownloadFiles()```. 

## Getting the information about a sample
The **information** about any **sample** may be retrived with something like:
```cpp
  //Select your dataset and load its information
  dm->LoadDataset("Sample_Name");
  
  //Now print some information
  cout << ">> Let's print some information..." << endl;
  //   + cross section
  cout << "   + X Section = " << dm->GetCrossSection() << endl;
  //   + Events in the sample
  cout << "   + N Events  = " << dm->GetEventsInTheSample() << endl;
  //   + Local Folder
  cout << "   + Local folder = " << dm->GetLocalFolder() << endl;
  //   + List of files
  cout << ">> Getting files..." << endl;
  std::vector<TString> files = dm->GetFiles();
  if (files.size() == 0) {
    cerr << "ERROR: Could not find files!" << endl;
    return;
  }

  //   + Dump all information to stdout
  cout << ">> Dumping..." << endl;
  dm->Dump();
``` 

## Dealing with Real Data Samples:
Finally, to get the files in a data sample the following syntax can be used:
```cpp
vector<TString> realdata= DatasetManager::GetRealDataFiles("DataName");
```

This call will return any file in the right place with a name matching the pattern
```[Tree_]DataName[_X[Y]].root``` where ```XY``` is a number.
For example the following names would match the example: ```DataName.root```, ```Tree_DataName.root```, 
```DataName_1.root```, ```DataName_01.root```, ```Tree_DataName_2.root```, ```Tree_DataName_12.root```, ...

# Full example

The file [UseExample.C](DatasetManager/UseExample.C) implements most of the code
above and can be used as an example of usage.

# Autorship
This class has been programmed by **I. González**. 

Refer to him if you find a bug or if you have proposals on how to extend the functionality.

Additional information about this class may be found in [(private) wiki][1].

[1]: http://www.hep.uniovi.es/wiki/index.php/DatasetManager