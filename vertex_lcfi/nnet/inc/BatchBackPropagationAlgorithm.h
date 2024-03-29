#ifndef BATCHBACKPROPAGATIONALGORITHM_H
#define BATCHBACKPROPAGATIONALGORITHM_H

#include "NeuralNetConfig.h"

#include "NeuralNet.h"
#include "NeuralNetDataSet.h"

#include <vector>

#ifdef __CINT__
#include "InputNormaliser.h"
#else
//namespace nnet added 15/08/06 by Mark Grimes (mark.grimes@bristol.ac.uk) for the LCFI vertex package
namespace nnet
{
class InputNormaliser;
}
#endif

//namespace nnet added 15/08/06 by Mark Grimes (mark.grimes@bristol.ac.uk) for the LCFI vertex package
namespace nnet
{

class
#ifndef __CINT__
NEURALNETDLL
#endif
BatchBackPropagationAlgorithm
{
public:
	BatchBackPropagationAlgorithm(NeuralNet &theNetwork,const double learningRate=0.5,const double momentumConstant=0.5);
	~BatchBackPropagationAlgorithm(void);
	BatchBackPropagationAlgorithm(const BatchBackPropagationAlgorithm&) = delete;
	BatchBackPropagationAlgorithm& operator=(const BatchBackPropagationAlgorithm&) = delete;
	void setLearningRate(const double newLearningRate)
	{ _learningRate = newLearningRate;}
	void setMomentumConstant(const double newMomentumConstant)
	{ _momentumConstant = newMomentumConstant;}
	void setMaxErrorIncrease(const double maxIncrease)
	{ _maxErrorInc = maxIncrease;}
	double train(const int numberOfEpochs,const NeuralNetDataSet &dataSet,
        const NeuralNet::InputNormalisationSelect normaliseTrainingData=NeuralNet::PassthroughNormalised);
    double train(const int numberOfEpochs,const NeuralNetDataSet &dataSet,const std::vector<InputNormaliser *> &inputNormalisers);
	void setProgressPrintoutFrequency(const int frequency) {_progressPrintoutFrequency = frequency;}
	void setEpochsToWaitBeforeRestore(const int epochs) {_epochsToWaitBeforeRestore = epochs;}
	std::vector<double> getTrainingErrorValuesPerEpoch() const {return _savedEpochErrorValues;}

protected:
    double trainWithDataSet(const int numberOfEpochs);
	std::vector<double> layerOutput(const int layer) const;
	void calculateLayerOutputs();
	void calculateDerivativeOutputs();
	void calculateErrorSignals();
	void calculateRunningGradientTotal();
	void calculateDeltaWeights();
	double error();
	double newEpoch();
	double processDataSet();

private:
	typedef std::vector<std::vector<double> > NetMatrix;

private:
	NeuralNet &_theNetwork;
	double _learningRate=0.0;
	double _momentumConstant=0.0;
	double _maxErrorInc=0.0;
	int _numberOfTrainingEvents=0;
	const std::vector<double> *_inputs=nullptr,*_target=nullptr;
	NetMatrix _neuronErrorSignals{};
	NetMatrix _neuronOutputs{};
	NetMatrix _neuronDerivativeOutputs{};
	NetMatrix _runningGradientTotal{};
	std::vector<double> _momentumWeights{};
	std::vector<double> _previousEpochWeights{};
	double _previousEpochError=0.0;
	double _runningEpochErrorTotal=0.0;
	const NeuralNetDataSet *_currentDataSet=nullptr;
	int _progressPrintoutFrequency=0;
	int _epochsToWaitBeforeRestore=0;
	std::vector<double> _savedEpochErrorValues{};
};

}//namespace nnet

#endif
