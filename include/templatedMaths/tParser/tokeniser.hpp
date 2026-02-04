#pragma once
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
//#include "tCmath.hpp"


/*****************************************\
!
!  This is a parsing library for parsing
!  mathematical expressions found in the
!  tCmath library and basic operators
!
\*****************************************/
struct Token{
  std::string type;
  std::string value;
  unsigned size;
};

/*****************************************\
!
!  A simple lexer that is used for
!  tokenization and lexographic analysis
!  expressions
!
\*****************************************/
void lex(std::string & InputStr, std::vector<Token> & tokens){
  std::stringstream InpSS;
  InpSS << InputStr;
  std::string parseVal;
  while(std::getline(InpSS, parseVal, ' ')){
    Token nTk;
    nTk.value = parseVal;
    tokens.push_back(nTk);
  };
};


/*****************************************\
!
!  A lexiographer that checks the type of
!  a value string containing a Var
!
\*****************************************/
void lexVarType(std::vector<Token> & VarTokens, std::vector<Token> & VarSizeTokens){
  for(unsigned I=0; I<VarTokens.size(); I++){
    unsigned Tensorflag=0, TRank=0;
    if(VarTokens[I].value.find("[") != std::string::npos) Tensorflag+=1;
    if(VarTokens[I].value.find("]") != std::string::npos) Tensorflag+=1;
    VarTokens[I].type = (Tensorflag==0) ? "scalar":((Tensorflag==2)?"tensor":"MALFTensor");
    if(Tensorflag!=0){
      TRank +=1;
      for(unsigned J=0; VarTokens[I].value[J] != 0; J++){
        if(VarTokens[I].value[J] == ',') TRank +=1;
      }
      for(unsigned J=0; J<VarSizeTokens.size(); J++){
        size_t Pos=VarTokens[I].value.find(VarSizeTokens[J].value);
//        if(Pos != std::string::npos) std::cout << std::setw(15) << Pos;
      }
    }
    VarTokens[I].size = TRank;
  }
};
