root -l /data/work/work_alice/JpsiReconstruction/AnalysisResults.root << eof
.L MRootIO.h
vector<string> vec_string
vector<TObject *> result = GetObjectRecursive(_file0,vec_string)
for (int i = 0; i < result.size(); i++) cout << result[i]->GetName() << endl;
cout << endl << endl;
for (int i = 0; i < vec_string.size(); i++) cout << vec_string[i] << endl;
eof
