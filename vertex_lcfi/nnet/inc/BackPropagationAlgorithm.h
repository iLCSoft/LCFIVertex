#ifndef BACKPROPAGATIONALGORITHM_H
#define BACKPROPAGATIONALGORITHM_H

#include "NeuralNetConfig.h"

#include "NeuralNet.h"

#include <vector>

//namespace nnet added 15/08/06 by Mark Grimes (mark.grimes@bristol.ac.uk) for the LCFI vertex package
namespace nnet
{

class
#ifndef __CINT__
NEURALNETDLL
#endif
BackPropagationAlgorithm
{
public:
	BackPropagationAlgorithm(NeuralNet &theNetwork,const double learningRate=0.5,const double momentumConstant=0.5);
	~BackPropagationAlgorithm(void);
	void setLearningRate(const double newLearningRate)
	{ _learningRate = newLearningRate;}
	void setMomentumConstant(const double newMomentumConstant)
	{ _momentumConstant = newMomentumConstant;}
	double train(const std::vector<double> &inputValues,const std::vector<double> &desiredOutput);
        BackPropagationAlgorithm(const BackPropagationAlgorithm&) = delete;
        BackPropagationAlgorithm& operator=(const BackPropagationAlgorithm&) = delete;
protected:
	std::vector<double> layerOutput(const int layer) const;
	void calculateLayerOutputs();
	void calculateDerivativeOutputs();
	void calculateErrorSignals();
	void calculateDeltaWeights();
	double error();

private:
	typedef std::vector<std::vector<double> > NetMatrix;

private:
	NeuralNet &_theNetwork;
	double _learningRate=0.0;
	const std::vector<double> *_inputs=nullptr,*_target=nullptr;
	NetMatrix _neuronErrorSignals{};
	NetMatrix _neuronOutputs{};
	NetMatrix _neuronDerivativeOutputs{};
	std::vector<double> _momentumWeights{};
	double _momentumConstant=0.0;
};

}//namespace nnet
#endif
