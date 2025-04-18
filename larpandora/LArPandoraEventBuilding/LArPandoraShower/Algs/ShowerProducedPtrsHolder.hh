//###################################################################
//### Name:        ShowerProducedPtrsHolder                       ###
//### Author:      Dominic Barker                                 ###
//### Date:        15.07.19                                       ###
//### Description: Class to holder the unique ptrs required to    ###
//###              produce data products. Used in                 ###
//###              LArPandoraModularShower and corresponding      ###
//###              tools.                                         ###
//###################################################################

#ifndef ShowerProducedPtrsHolder_HH
#define ShowerProducedPtrsHolder_HH

#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Algs/ShowerElementHolder.hh"

//Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/demangle.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace reco::shower {
  class ShowerUniqueProduerPtrBase;
  template <class T>
  class ShowerUniqueProductPtr;
  template <class T>
  class ShowerUniqueAssnPtr;
  class ShowerPtrMakerBase;
  template <class T>
  class ShowerPtrMaker;
  class ShowerProducedPtrsHolder;
}

#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <utility> // std::move()
#include <vector>

//Template to check if its an assn
template <class N>
struct is_assn {
  static const int value = 0;
};

template <class L, class R, class D>
struct is_assn<art::Assns<L, R, D>> {
  static const int value = 1;
};

//Template to carry the type becuase functions can't, annoyingly, be partially specialised
template <typename T>
struct type {};

//Base class for the class that holds the ptr.
class reco::shower::ShowerUniqueProduerPtrBase {

public:
  virtual ~ShowerUniqueProduerPtrBase() noexcept = default;

  virtual void reset() = 0;

  virtual void AddDataProduct(const reco::shower::ShowerElementHolder& selement_holder,
                              const std::string& Name) = 0;

  virtual void MoveToEvent(art::Event& evt) = 0;

  virtual std::string GetType() const = 0;
  virtual std::string GetInstanceName() const = 0;

  virtual int GetVectorPtrSize() const { return -1; }
};

//Class that holds a unique ptr for the product. This is what is stored in the map. The product is put into
//the event as a vector so the this holder maintains this unique ptr and the actions to manipulate it.
template <class T>
class reco::shower::ShowerUniqueProductPtr<std::vector<T>>
  : public reco::shower::ShowerUniqueProduerPtrBase {

public:
  ShowerUniqueProductPtr(const std::string& Instancename)
  {
    ptr = 1;
    showeruniqueptr = std::make_unique<std::vector<T>>();
    InstanceName = Instancename;
  }

  //Get the unique ptr for the data product.
  std::unique_ptr<T>& GetPtr()
  {
    if (ptr) { return showeruniqueptr; }
    else {
      throw cet::exception("ShowerUniqueProduerPtr") << "Element does not exist" << std::endl;
    }
  }

  void reset() override { showeruniqueptr.reset(new std::vector<T>()); }

  //Add a data product on to the vector that will be added to the event.
  void AddDataProduct(const reco::shower::ShowerElementHolder& selement_holder,
                      const std::string& Name) override
  {
    T product;
    if (!selement_holder.CheckElement(Name)) {
      mf::LogError("ShowerProducedPtrsHolder")
        << "Trying to add data product: " << Name
        << ". This element does not exist in the element holder" << std::endl;
      return;
    }

    int err = selement_holder.GetElement(Name, product);
    if (err) {
      mf::LogError("ShowerProducedPtrsHolder")
        << "Trying to add data product: " << Name
        << ". This element does not exist in the element holder" << std::endl;
      return;
    }
    showeruniqueptr->push_back(product);
    return;
  }

  //Final thing to do move to the event.
  void MoveToEvent(art::Event& evt) override { evt.put(std::move(showeruniqueptr), InstanceName); }

  //This returns the type of the undrlying element. It should be of the form std::unique_ptr<std::vector<T> >
  std::string GetType() const override
  {
    return cet::demangle_symbol(typeid(showeruniqueptr.get()).name());
  }

  //A user can set the instances name like an other producer module product. Get it using this.
  std::string GetInstanceName() const override { return InstanceName; }

  //Get the size of the std::vector. Usful for making assns.
  int GetVectorPtrSize() const override { return showeruniqueptr->size(); }

private:
  //Element itself that is put into the art event.
  std::unique_ptr<std::vector<T>> showeruniqueptr;

  //bool to see if the element is set.
  bool ptr;

  //Name when saved into the the event default is ""
  std::string InstanceName;
};

//Class that holds a unique ptr for the association. This is what is stored in the map. The association is put into
//the event as a vector so the this holder maintains this unique ptr and the actions to manipulate it.
//I guess if the product if a product is unique to the event then this holder will deal with it.
//I have tried to be smart and I don't think it not being an association is a problem as
//long as the user is smart.
template <class T>
class reco::shower::ShowerUniqueAssnPtr : public reco::shower::ShowerUniqueProduerPtrBase {

public:
  ShowerUniqueAssnPtr(const std::string& Instancename)
  {
    ptr = 1;
    showeruniqueptr = std::make_unique<T>();
    InstanceName = Instancename;
  }

  //Get the ptr to the association.
  std::unique_ptr<T>& GetPtr()
  {
    if (ptr) { return showeruniqueptr; }
    else {
      throw cet::exception("ShowerUniqueAssnPtr") << "Element does not exist" << std::endl;
    }
  }

  void reset() override { showeruniqueptr.reset(new T()); }

  //place the association to the event.
  void MoveToEvent(art::Event& evt) override { evt.put(std::move(showeruniqueptr), InstanceName); }

  //Not need but the compiler complains if its not here.
  void AddDataProduct(const reco::shower::ShowerElementHolder& selement_holder,
                      const std::string& Name) override
  {
    throw cet::exception("ShowerUniqueAssnPtr")
      << "The creator of this code has failed you. Please contact Dominic Bakrer" << std::endl;
  }

  //Get the type in form as a string. It should be the form of std::unique_ptr<art::Assn<T A, T1 B> >
  std::string GetType() const override
  {
    return cet::demangle_symbol(typeid(showeruniqueptr.get()).name());
  }

  //Get the Instance name that product will be saved as in the art::Event.
  std::string GetInstanceName() const override { return InstanceName; }

private:
  //Actual Element the assn. This is put into the event.
  std::unique_ptr<T> showeruniqueptr;

  //bool to see if the element is filled.
  bool ptr;

  //Name when saved into the the event default is ""
  std::string InstanceName;
};

//Base class to hold the pointer makers. This interacts the with the module a little differently
//as the ptr makers do not work within the tools. The holds and the ptrmaker and provides set and
//get functions so that the usr can access art::Ptrs for the products of the module and associate
//them to other things.
class reco::shower::ShowerPtrMakerBase {

public:
  virtual ~ShowerPtrMakerBase() noexcept = default;

  virtual bool CheckPtrMaker() const = 0;

  virtual void SetPtrMaker(art::Event& evt) = 0;

  virtual void Reset() = 0;
};

//Derived class - See above
template <class T>
class reco::shower::ShowerPtrMaker : public ShowerPtrMakerBase {

public:
  ShowerPtrMaker(const std::string& Instancename)
  {
    ptrmaker = nullptr;
    ptr = 0;
    InstanceName = Instancename;
  }

  //Check the ptr maker is ready to be used.
  bool CheckPtrMaker() const override
  {
    if (ptr) { return true; }
    return false;
  }

  //Return the ptr maker. Probably never needed.
  art::PtrMaker<T>& GetPtrMaker()
  {
    if (ptr) {
      if (ptrmaker == nullptr) {
        throw cet::exception("ShowerPtrMaker") << "Ptr maker ptr is null" << std::endl;
      }
      return *ptrmaker;
    }
    throw cet::exception("ShowerPtrMaker")
      << "Trying to get a  ptrmaker that does not exists" << std::endl;
  }

  //Return the art ptr that the module produces corresponding the index iter
  art::Ptr<T> GetArtPtr(int iter) const
  {
    if (ptr) {
      if (ptrmaker == nullptr) {
        throw cet::exception("ShowerPtrMaker") << "Ptr maker ptr is null" << std::endl;
      }
      return (*ptrmaker)(iter);
    }
    throw cet::exception("ShowerPtrMaker")
      << "Trying to get a  ptrmaker that does not exists" << std::endl;
  }

  //Set the ptr maker this is reset at the start of the event.
  void SetPtrMaker(art::Event& evt) override
  {
    ptrmaker.reset(new art::PtrMaker<T>(evt, InstanceName));
    ptr = 1;
  }

  void Reset() override
  {
    if (!ptr) {
      throw cet::exception("ShowerPtrMaker") << "Trying to reset ptr but it has not been set in "
                                                "the first place. Please contatc Dom Barker"
                                             << std::endl;
    }
    ptrmaker.reset(nullptr);
    ptr = 0;
  }

private:
  //The ptr maker itself. Used to make art::Ptrs to make assns.
  std::unique_ptr<art::PtrMaker<T>> ptrmaker;

  //The name of the data product which will be saved in the event. The ptr maker requires this.
  std::string InstanceName;

  //bool to tell if the ptr maker is ready to use or not.
  int ptr;
};

//Class that holds all the unique ptrs and the ptr makers. It is what the tools and module use
//to access the above class elements. The end user case will see the user not interact with this
// directly.
class reco::shower::ShowerProducedPtrsHolder {

public:
  //Initialise the a unique ptr in the map. This will be added to the event.
  template <class T>
  int SetShowerUniqueProduerPtr(type<T>, const std::string& Name, const std::string& Instance = "")
  {

    //Add to the assns
    if (showerassnPtrs.find(Name) != showerassnPtrs.end()) {
      mf::LogWarning("ShowerProducedPtrsHolder")
        << "Trying to set Element: " << Name << ". This element has already been set. Please check."
        << std::endl;
      return 1;
    }

    //Check the same type has not already been set.
    if (!CheckForMultipleTypes(type<T>(), Name, Instance)) {
      throw cet::exception("ShowerProducedPtrsHolder")
        << "Trying to set multiple objects with same type with no instance name or same instance "
           "name"
        << std::endl;
    }

    showerassnPtrs[Name] = std::make_unique<ShowerUniqueAssnPtr<T>>(Instance);
    return 0;
  }

  //Set the unique ptr. The unique ptr will be filled into the event.
  template <class T>
  int SetShowerUniqueProduerPtr(type<std::vector<T>>,
                                const std::string& Name,
                                const std::string& Instance = "")
  {

    //Then add the products
    if (showerproductPtrs.find(Name) != showerproductPtrs.end()) {
      mf::LogWarning("ShowerProducedPtrsHolder")
        << "Trying to set Element: " << Name << ". This element has already been set. Please check."
        << std::endl;
      return 1;
    }

    //Check the same type has not already been set.
    if (!CheckForMultipleTypes(type<std::vector<T>>(), Name, Instance)) {
      throw cet::exception("ShowerProducedPtrsHolder")
        << "Trying to set multiple objects with same type with no instance name or same instance "
           "name"
        << std::endl;
    }

    if (showerPtrMakers.find(Name) != showerPtrMakers.end()) {
      throw cet::exception("ShowerProducedPtrsHolder")
        << "PtrMaker already exist. It should not be set again" << std::endl;
    }
    showerPtrMakers[Name] = std::make_unique<ShowerPtrMaker<T>>(Instance);
    showerproductPtrs[Name] = std::make_unique<ShowerUniqueProductPtr<std::vector<T>>>(Instance);
    return 0;
  }

  //Checks if the ptr exists
  bool CheckUniqueProduerPtr(const std::string& Name) const
  {
    if (showerproductPtrs.find(Name) != showerproductPtrs.end()) { return true; }
    if (showerassnPtrs.find(Name) != showerassnPtrs.end()) { return true; }
    return false;
  }

  //Reset the ptrs;
  void reset()
  {
    for (auto const& showerptr : showerproductPtrs) {
      (showerptr.second)->reset();
    }
    for (auto const& showerptr : showerassnPtrs) {
      (showerptr.second)->reset();
    }
  }

  //Add any data products that are produced by the module to the unique ptr it corresponds to
  //This is done by matching strings in the element holder and the ptr holder. Hence these
  //must match. This is a global command done in the module.
  void AddDataProducts(const reco::shower::ShowerElementHolder& selement_holder)
  {
    for (auto const& showerproductPtr : showerproductPtrs) {
      (showerproductPtr.second)->AddDataProduct(selement_holder, showerproductPtr.first);
    }
  }

  //Global command to move all products into the event. This is done in the module.
  void MoveAllToEvent(art::Event& evt)
  {
    for (auto const& showerproductPtr : showerproductPtrs) {
      (showerproductPtr.second)->MoveToEvent(evt);
    }
    for (auto const& showerassnPtr : showerassnPtrs) {
      (showerassnPtr.second)->MoveToEvent(evt);
    }
  }

  bool CheckAllProducedElements(reco::shower::ShowerElementHolder& selement_holder) const
  {
    bool checked = true;
    for (auto const& showerproductPtr : showerproductPtrs) {
      if (showerproductPtr.first == "shower") { continue; }
      checked = checked && selement_holder.CheckElement(showerproductPtr.first);
    }
    return checked;
  }

  //This returns the unique ptr. This is a legacy code.
  template <class T>
  T& GetPtr(const std::string& Name)
  {
    auto const showerproductPtrsIt = showerproductPtrs.find(Name);
    if (showerproductPtrsIt != showerproductPtrs.end()) {
      reco::shower::ShowerUniqueProductPtr<T>* prod =
        dynamic_cast<reco::shower::ShowerUniqueProductPtr<T>*>(showerproductPtrsIt->second.get());
      return prod->GetPtr();
    }

    auto const showerassnPtrsIt = showerassnPtrs.find(Name);
    if (showerassnPtrsIt != showerassnPtrs.end()) {
      reco::shower::ShowerUniqueAssnPtr<T>* assn =
        dynamic_cast<reco::shower::ShowerUniqueAssnPtr<T>*>(showerassnPtrsIt->second.get());
      return assn->GetPtr();
    }

    throw cet::exception("ShowerProducedPtrsHolder")
      << "Trying to get Ptr for: " << Name << " but Element does not exist" << std::endl;
  }

  //Wrapper so that the use the addSingle command for the association. Add A and B to the association just
  //as if add single add.
  template <class T, class A, class B>
  void AddSingle(A& a, B& b, const std::string& Name)
  {
    auto const showerassnPtrsIt = showerassnPtrs.find(Name);
    if (showerassnPtrsIt == showerassnPtrs.end()) {
      throw cet::exception("ShowerProducedPtrsHolder")
        << "Trying to get the association: " << Name << "Element does not exist" << std::endl;
    }
    if (!is_assn<T>::value) {
      throw cet::exception("ShowerProducedPtrsHolder")
        << "Element type  is not an assoication please only use this for assocations" << std::endl;
    }
    reco::shower::ShowerUniqueAssnPtr<T>* assnptr =
      dynamic_cast<reco::shower::ShowerUniqueAssnPtr<T>*>(showerassnPtrsIt->second.get());
    if (assnptr == nullptr) {
      throw cet::exception("ShowerProducedPtrsHolder")
        << "Failed to cast back. Maybe you got the type wrong or you are accidently accessing a "
           "differently named product"
        << std::endl;
    }

    T* assn = dynamic_cast<T*>(assnptr->GetPtr().get());
    if (assn == nullptr) {
      throw cet::exception("ShowerProducedPtrsHolder")
        << "Something went wrong trying to cast tothe assn. Maybe the name: " << Name
        << " exists but its not an assn" << std::endl;
    }

    assn->addSingle(a, b);
    return;
  }

  //Initialise the ptr makers. This is done at the the start of the module.
  void SetPtrMakers(art::Event& evt)
  {
    for (auto const& showerPtrMaker : showerPtrMakers) {
      if (showerPtrMakers.find(showerPtrMaker.first) == showerPtrMakers.end()) {
        throw cet::exception("ShowerProducedPtrsHolder")
          << "PtrMaker was empty. This is concerning" << std::endl;
      }
      showerPtrMakers[showerPtrMaker.first]->SetPtrMaker(evt);
    }
  }

  //Wrapper to access a particle PtrMaker. This is legacy as is not used.
  template <class T>
  art::PtrMaker<T>& GetPtrMaker(const std::string& Name)
  {
    auto const showerPtrMakersIt = showerPtrMakers.find(Name);
    if (showerPtrMakersIt == showerPtrMakers.end()) {
      throw cet::exception("ShowerProducedPtrsHolder") << "PtrMaker does not exist" << std::endl;
    }
    else {
      if (!showerPtrMakersIt->second->CheckPtrMaker()) {
        throw cet::exception("ShowerProducedPtrsHolder") << "PtrMaker is not set" << std::endl;
      }
      reco::shower::ShowerPtrMaker<T>* ptrmaker =
        dynamic_cast<reco::shower::ShowerPtrMaker<T>*>(showerassnPtrs[Name].get());
      return ptrmaker->GetPtrMaker();
    }
  }

  //Wrapper to return to the the user the art ptr corresponding the index iter
  template <class T>
  art::Ptr<T> GetArtPtr(const std::string& Name, const int& iter) const
  {
    auto const showerPtrMakersIt = showerPtrMakers.find(Name);
    if (showerPtrMakersIt == showerPtrMakers.end()) {
      throw cet::exception("ShowerProducedPtrsHolder")
        << "PtrMaker does not exist for " << Name << " Did you initialise this? " << std::endl;
    }
    else {
      if (!showerPtrMakersIt->second->CheckPtrMaker()) {
        throw cet::exception("ShowerProducedPtrsHolder")
          << "PtrMaker is not set. This is an issue for the devlopment team me. Contact Dom Barker"
          << std::endl;
      }
      reco::shower::ShowerPtrMaker<T>* ptrmaker =
        dynamic_cast<reco::shower::ShowerPtrMaker<T>*>(showerPtrMakersIt->second.get());
      if (ptrmaker == nullptr) {
        throw cet::exception("ShowerProducedPtrsHolder")
          << "Failed to cast back. Maybe you got the type wrong or you are accidently accessing a "
             "differently named product"
          << std::endl;
      }
      return ptrmaker->GetArtPtr(iter);
    }
  }

  //Legacy not used.
  void ResetPtrMakers()
  {
    for (auto const& showerPtrMaker : showerPtrMakers) {
      (showerPtrMaker.second)->Reset();
    }
  }

  //Return the size of the std::vector of the data object with the unique name string.
  int GetVectorPtrSize(const std::string& Name) const
  {
    auto const showerproductPtrsIt = showerproductPtrs.find(Name);
    if (showerproductPtrsIt != showerproductPtrs.end()) {
      return showerproductPtrsIt->second->GetVectorPtrSize();
    }
    throw cet::exception("ShowerProducedPtrsHolder")
      << "Product: " << Name << " has not been set in the producers map" << std::endl;
  }

  //Print the type, the element name and the instance name
  void PrintPtr(const std::string& Name) const
  {
    auto const showerproductPtrsIt = showerproductPtrs.find(Name);
    if (showerproductPtrsIt != showerproductPtrs.end()) {
      const std::string Type = showerproductPtrsIt->second->GetType();
      const std::string InstanceName = showerproductPtrsIt->second->GetInstanceName();
      std::cout << "Element Name: " << Name << " Instance Name: " << InstanceName
                << " Type: " << Type << std::endl;
      return;
    }
    auto const showerassnPtrsIt = showerassnPtrs.find(Name);
    if (showerassnPtrsIt != showerassnPtrs.end()) {
      const std::string Type = showerassnPtrsIt->second->GetType();
      const std::string InstanceName = showerassnPtrsIt->second->GetInstanceName();
      std::cout << "Element Name: " << Name << " Instance Name: " << InstanceName
                << " Type: " << Type << std::endl;
      return;
    }
    mf::LogError("ShowerProducedPtrsHolder")
      << "Trying to print Element: " << Name
      << ". This element does not exist in the element holder" << std::endl;
    return;
  }

  //This function will print out all the pointers and there types for the user to check.
  void PrintPtrs() const
  {

    unsigned int maxname = 0;
    for (auto const& showerprodPtr : showerproductPtrs) {
      if (showerprodPtr.first.size() > maxname) { maxname = showerprodPtr.first.size(); }
    }
    for (auto const& showerassnPtr : showerassnPtrs) {
      if (showerassnPtr.first.size() > maxname) { maxname = showerassnPtr.first.size(); }
    }

    std::map<std::string, std::pair<std::string, std::string>> Type_showerprodPtrs;
    std::map<std::string, std::pair<std::string, std::string>> Type_showerassnPtrs;
    for (auto const& showerprodPtr : showerproductPtrs) {
      const std::string Type = (showerprodPtr.second)->GetType();
      const std::string InstanceName = (showerprodPtr.second)->GetInstanceName();
      Type_showerprodPtrs[showerprodPtr.first] = std::make_pair(InstanceName, Type);
    }
    for (auto const& showerassnPtr : showerassnPtrs) {
      const std::string Type = (showerassnPtr.second)->GetType();
      const std::string InstanceName = (showerassnPtr.second)->GetInstanceName();
      Type_showerassnPtrs[showerassnPtr.first] = std::make_pair(InstanceName, Type);
    }

    unsigned int maxtype = 0;
    unsigned int maxinstname = 0;
    for (auto const& Type_showerprodPtr : Type_showerprodPtrs) {
      if (Type_showerprodPtr.second.second.size() > maxtype) {
        maxtype = Type_showerprodPtr.second.second.size();
      }
      if (Type_showerprodPtr.second.first.size() > maxinstname) {
        maxinstname = Type_showerprodPtr.second.first.size();
      }
    }
    for (auto const& Type_showerassnPtr : Type_showerassnPtrs) {
      if (Type_showerassnPtr.second.second.size() > maxtype) {
        maxtype = Type_showerassnPtr.second.second.size();
      }
      if (Type_showerassnPtr.second.first.size() > maxinstname) {
        maxinstname = Type_showerassnPtr.second.first.size();
      }
    }

    unsigned int n = maxname + maxtype + maxinstname + 51;
    std::cout << std::left << std::setfill('*') << std::setw(n - 1) << "**" << std::endl;
    std::cout << "Unique Ptrs that are added to the event" << std::endl;
    std::cout << std::left << std::setfill('*') << std::setw(n - 1) << "**" << std::endl;
    for (auto const& Type_showerprodPtr : Type_showerprodPtrs) {
      std::cout << std::left << std::setfill(' ') << std::setw(21)
                << "* Data Product Name: " << std::setw(maxname) << Type_showerprodPtr.first;
      std::cout << std::left << std::setfill(' ') << " * Instance Name: " << std::setw(maxinstname)
                << Type_showerprodPtr.second.first;
      std::cout << std::left << std::setfill(' ') << " * Type: " << std::setw(maxtype)
                << Type_showerprodPtr.second.second << " *" << std::endl;
    }
    for (auto const& Type_showerassnPtr : Type_showerassnPtrs) {
      std::cout << std::left << std::setfill(' ') << std::setw(maxname) << std::setw(21)
                << "* Association Name: " << std::setw(maxname) << Type_showerassnPtr.first;
      std::cout << std::left << std::setfill(' ') << " * Instance Name: " << std::setw(maxinstname)
                << Type_showerassnPtr.second.first;
      std::cout << std::left << std::setfill(' ') << " * Type: " << std::setw(maxtype)
                << Type_showerassnPtr.second.second << " *" << std::endl;
    }
    std::cout << std::left << std::setfill('*') << std::setw(n - 1) << "**" << std::endl;
    std::cout << std::setfill(' ');
    std::cout << std::setw(0);
    return;
  }

private:
  //Function to check that a data product does not already exist with the same instance name
  template <class T>
  bool CheckForMultipleTypes(type<T>, const std::string& Name, const std::string& Instance) const
  {

    //Check the a product of the same does not exist without a different instance name
    for (auto const& assn : showerassnPtrs) {
      reco::shower::ShowerUniqueAssnPtr<T>* assnptr =
        dynamic_cast<reco::shower::ShowerUniqueAssnPtr<T>*>(assn.second.get());
      if (assnptr != nullptr) {
        if (assnptr->GetInstanceName() == Instance) { return false; }
      }
    }
    return true;
  }

  //Function to check that a data product does not already exist with the same instance name
  template <class T>
  bool CheckForMultipleTypes(type<std::vector<T>>,
                             const std::string& Name,
                             const std::string& Instance) const
  {

    //Check the a product of the same does not exist without a different instance name
    for (auto const& product : showerproductPtrs) {
      reco::shower::ShowerUniqueProductPtr<std::vector<T>>* prod =
        dynamic_cast<reco::shower::ShowerUniqueProductPtr<std::vector<T>>*>(product.second.get());
      if (prod != nullptr) {
        if (prod->GetInstanceName() == Instance) { return false; }
      }
    }
    return true;
  }

  //Holder of the data objects of type std::vector<T> that will be saved in the events.
  std::map<std::string, std::unique_ptr<reco::shower::ShowerUniqueProduerPtrBase>>
    showerproductPtrs;

  //Holder of the data objects of type T that will be saved into the events. I think these will only be assns.
  std::map<std::string, std::unique_ptr<reco::shower::ShowerUniqueProduerPtrBase>> showerassnPtrs;

  //Holder of the ptrMakers whcih make the art::Ptrs of the data objects that lie in showerproductPtrs.
  std::map<std::string, std::unique_ptr<reco::shower::ShowerPtrMakerBase>> showerPtrMakers;
};

#endif
