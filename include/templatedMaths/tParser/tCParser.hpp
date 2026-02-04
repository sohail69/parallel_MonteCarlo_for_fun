#pragma once
#include <map>
#include <vector>
#include <string>
#include <functional>
//#include "tCmath.hpp"
#include "tokeniser.hpp"

/*****************************************\
!
!  This is the parsing tree basic Node
!  NodeType=="Operator"
!  NodeType=="Function"
!  NodeType=="Var"
!
\*****************************************/
struct BTreeNode{
  std::string Data;
  std::string NodeType;
  BTreeNode* LNode=NULL;
  BTreeNode* RNode=NULL;
};


/*****************************************\
!
!  A simple parser that is used for
!  parsing tCmath expressions with
!  tensors
!
\*****************************************/
template<typename Number>
std::function<Number(std::vector<Number> data)> tensorParse(std::string &Iters
                                                          , std::string &varSizes
                                                          , std::string &Vars
                                                          , std::string &expr)
{
  std::vector<Token> IterTKNs, varSizeTKNs, varTKNs;
  lex(Iters, IterTKNs);

  lex(varSizes, varSizeTKNs);

  lex(Vars, varTKNs);
  lexVarType(varTKNs, varSizeTKNs);

  auto lmbda2 = [=](std::vector<Number> data){
    Number a(0.00);
    for(unsigned I=0; I<10; I++){
      a = a + Number(1.00);
    }
    return a;
  };
  return lmbda2;
};
