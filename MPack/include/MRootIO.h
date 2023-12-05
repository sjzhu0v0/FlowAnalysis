#ifndef MROOTIO_H
#define MROOTIO_H

#include "MSystem.h"

vector<TObject *> GetObjectRecursive(TObject *folder,
                                     vector<string> &vec_string) {
  vector<TObject *> vec_obj;

  if (!folder)
    return vec_obj;

  if (!folder->IsFolder()) {
    vec_obj.push_back(folder);
    vec_string.push_back(folder->GetName());
    return vec_obj;
  }

  // TCollection TDirecotry
  if (folder->IsA()->InheritsFrom(TDirectory::Class())) {
    TDirectory *dir = static_cast<TDirectory *>(folder);
    TIter next(dir->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next())) {
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom(TObject::Class()))
        continue;
      TObject *obj = key->ReadObj();
      vector<string> vec_string_tmp;
      vector<TObject *> vec_obj_tmp = GetObjectRecursive(obj,vec_string_tmp);
      for (auto &str : vec_string_tmp) {
        str = (string)folder->GetName() + "/" + str;
      }
      vec_obj.insert(vec_obj.end(), vec_obj_tmp.begin(), vec_obj_tmp.end());
      vec_string.insert(vec_string.end(), vec_string_tmp.begin(),
                        vec_string_tmp.end());
    }
  } else if (folder->IsA()->InheritsFrom(TList::Class())) {
    TList *list = static_cast<TList *>(folder);
    TIter next(list);
    TObject *obj;
    while ((obj = next())) {
      vector<string> vec_string_tmp;
      vector<TObject *> vec_obj_tmp = GetObjectRecursive(obj,vec_string_tmp);
      for (auto &str : vec_string_tmp) {
        str = (string)folder->GetName() + "/" + str;
      }
      vec_obj.insert(vec_obj.end(), vec_obj_tmp.begin(), vec_obj_tmp.end());
      vec_string.insert(vec_string.end(), vec_string_tmp.begin(),
                        vec_string_tmp.end());
    }
  } else {
    cout << "Unknown class: " << folder->ClassName() << endl;
  }

  return vec_obj;
}

template <typename T>
vector<T *> GetObjectVector(TDirectory *dir = gDirectory) {
  vector<T *> list;
  TIter next(dir->GetListOfKeys());
  TKey *key;
  while ((key = (TKey *)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom(T::Class()))
      continue;
    T *obj = (T *)key->ReadObj();
    list.push_back(obj);
  }
  return list;
}

template <typename T>
vector<T *> GetObjectVector(vector<TObject*> vec_obj, vector<string> vec_string_input, vector<string> &vec_string_output, string str_tag = "") {
  vector<T *> list;
  for (int i = 0; i < vec_string_input.size(); i++) {
    if (vec_obj[i]->InheritsFrom(T::Class())) {
      if (strcmp(str_tag.c_str(), "") == 0) {
        list.push_back((T *)vec_obj[i]);
        vec_string_output.push_back(vec_string_input[i]);
      } else {
        if (vec_string_input[i].find(str_tag) != string::npos) {
          list.push_back((T *)vec_obj[i]);
          vec_string_output.push_back(vec_string_input[i]);
        }
      }
    }
  }
  return list;
}

template <typename T>
vector<T *> GetObjectVector(vector<TObject*> vec_obj, vector<string> vec_string_input, vector<string> vec_string_output, TString str_tag = "") {
  return GetObjectVector<T>(vec_obj, vec_string_input, vec_string_output, (string)str_tag);
}

std::vector<std::vector<void *>> GetDataAddressesFromTree(TTree *tree) {
  std::vector<std::vector<void *>> dataAddresses;
  for (auto branch : *tree->GetListOfBranches()) {
    std::vector<void *> branchAddresses;
    for (auto leaf : *static_cast<TBranch *>(branch)->GetListOfLeaves()) {
      branchAddresses.push_back(static_cast<TLeaf *>(leaf)->GetValuePointer());
    }
    dataAddresses.push_back(branchAddresses);
  }
  return dataAddresses;
}

template <typename T> T GetValueFromAddress(void *address) {
  return *(static_cast<T *>(address));
}

#endif