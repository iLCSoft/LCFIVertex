#ifndef INPUTIMPORTANCE_H
#define INPUTIMPORTANCE_H

#include "NeuralNetConfig.h"
#include "NeuralNet.h"
#include "NeuralNetDataSet.h"

#include <vector>

class
#ifndef __CINT__
NEURALNETDLL
#endif
InputImportance
{
public:
	std::vector<double> operator()(const NeuralNet &theNet,const NeuralNetDataSet &theData) const;
};
#endif
