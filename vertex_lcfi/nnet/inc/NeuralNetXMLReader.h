#ifndef NEURALNETXMLREADER_H
#define NEURALNETXMLREADER_H
#ifndef NEURALNETNOXMLREADER
#include "NeuralNetConfig.h"

#if defined(_WIN32)&&!defined(USEXERCES)&&!defined(USETINYXML)
# import <msxml4.dll>
#elif defined(USEXERCES)&&!defined(USETINYXML)
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#else
#ifndef TIXML_USE_STL // TINYXML using STL classes
#define TIXML_USE_STL
#endif
//Change: 15/08/06 by Mark Grimes (mark.grimes@bristol.ac.uk)
//Added 'marlin/' here for the LCFI vertex package to reflect the Marlin directory structure
#include "marlin/tinyxml.h"
#endif

#include <string>
#include <vector>

//namespace nnet added 15/08/06 by Mark Grimes (mark.grimes@bristol.ac.uk) for the LCFI vertex package
namespace nnet
{

class NEURALNETDLL NeuralNetXMLReader
{
public:
	NeuralNetXMLReader(const std::string &xmlFileName);
	NeuralNetXMLReader(const char *xmlFileName);
	~NeuralNetXMLReader();
	NeuralNetXMLReader(const nnet::NeuralNetXMLReader&) = delete;
	NeuralNetXMLReader& operator=(const nnet::NeuralNetXMLReader&) = delete;
	void networkTopology(int &numberOfInputs,int &numberOfLayers);
    void getNormalisationData(std::vector<std::string> &normaliserTypes,std::vector<std::vector<double> > &constructionData,
        		              std::vector<double> &targetoffsets,std::vector<double> &targetranges);
	int nextLayer();
	int nextNeuron(std::string &neuronType,std::vector<double> &constructionData);

private:
	void init();
	void clean();
	void getLayerCollection();
	void getNeuronCollection();

private:
	std::string _xmlFileName{};
	bool _ok=false;
	int _currentLayer=0;
	int _currentNeuron=0;
#if defined(_WIN32)&&!defined(USEXERCES)&&!defined(USETINYXML)
	MSXML2::IXMLDOMDocumentPtr _pXMLDoc{};
	MSXML2::IXMLDOMNodeListPtr _pLayers{};
	MSXML2::IXMLDOMNodeListPtr _pNeurons{};
#elif defined(USEXERCES)&&!defined(USETINYXML)
	XERCES_CPP_NAMESPACE_QUALIFIER XercesDOMParser *_pDOMParser=nullptr;
	XERCES_CPP_NAMESPACE_QUALIFIER ErrorHandler *_errorHandler=nullptr;
	XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument *_pDoc=nullptr;
	XERCES_CPP_NAMESPACE_QUALIFIER DOMNodeList *_pLayers=nullptr;
	XERCES_CPP_NAMESPACE_QUALIFIER DOMNodeList *_pNeurons=nullptr;
#else
    TiXmlDocument *_pXMLDoc=nullptr;
    TiXmlElement *_pDocRoot=nullptr;
    std::vector<TiXmlNode *> *_pLayers=nullptr;
    std::vector<TiXmlNode *> *_pNeurons=nullptr;
#endif
};

}//namespace nnet

#endif
#endif
