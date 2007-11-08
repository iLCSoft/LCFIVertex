#Basic example python script for Jas3
#Makes plots stacked plots of the FlavourTagInputs
#First open LCFIPlot.aida and this file, and then "F2" to run.
#Uncomment the lines at the bottom to make plots of different variables - more than one may be made at once
#
#Ignore the error: trying to follow the advice made the problem worse
#"There might be a synchronization problem when the plotter is being shown and it is written to file.
# Please refer to the JAIDA release notes for more details: http://java.freehep.org/jaida.
# To avoid this problem don't invoke the show() method on the Plotter when writing to file."
#
#victoria.martin@ed.ac.uk

from hep.aida import *
from java.util import Random
from java.lang import Boolean

true = Boolean("true");
false = Boolean("false");
    
def FlavourTagInputsOverlay(hname, htitle):
    "plot bjet, cjet, udsjet histograms overlayed"

    factory = IAnalysisFactory.create();
    tree = factory.createTreeFactory().create();
    hf = factory.createHistogramFactory(tree);
        
    plotterFactory = factory.createPlotterFactory();
    plotter = plotterFactory.create("Flavour Tag Inputs");
    
    bjet = aidaMasterTree.find("/LCFIPlot.aida/MyLCFIAIDAPlotProcessor/FlavourTagInputs/bJets/"+hname);
    cjet = aidaMasterTree.find("/LCFIPlot.aida/MyLCFIAIDAPlotProcessor/FlavourTagInputs/cJets/"+hname);
    ljet = aidaMasterTree.find("/LCFIPlot.aida/MyLCFIAIDAPlotProcessor/FlavourTagInputs/udsJets/"+hname);
 
    tree.mkdir("/LCFI Flavour Tag Inputs: " + hname);
    tree.cd("/LCFI Flavour Tag Inputs: " + hname);
  
    bjet.scale((1./bjet.sumAllBinHeights()));
    cjet.scale((1./cjet.sumAllBinHeights()));
    ljet.scale((1./ljet.sumAllBinHeights()));

    lSum = hf.createCopy("Light Jets",bjet);
    clSum = hf.add("C Jets",cjet,lSum);
    bclSum = hf.add("B Jets",bjet,clSum);
 
    normFactor = bclSum.sumAllBinHeights();
  
    bclSum.scale(1./normFactor);
    lSum.scale(1./normFactor);
    clSum.scale(1./normFactor);
        
    bDataStyle = plotterFactory.createPlotterStyle();
    cDataStyle = plotterFactory.createPlotterStyle();
    lDataStyle = plotterFactory.createPlotterStyle();
    
    bDataStyle.dataStyle().lineStyle().setParameter("color","red")
    cDataStyle.dataStyle().lineStyle().setParameter("color","green")
    lDataStyle.dataStyle().lineStyle().setParameter("color","blue")   
    bDataStyle.dataStyle().lineStyle().setParameter("thickness","2")
    cDataStyle.dataStyle().lineStyle().setParameter("thickness","2")
    lDataStyle.dataStyle().lineStyle().setParameter("thickness","2")
     
    bDataStyle.dataStyle().lineStyle().setVisible(true);
    cDataStyle.dataStyle().lineStyle().setVisible(true);
    lDataStyle.dataStyle().lineStyle().setVisible(true);
    
    bDataStyle.dataStyle().markerStyle().setVisible(false);
    cDataStyle.dataStyle().markerStyle().setVisible(false);
    lDataStyle.dataStyle().markerStyle().setVisible(false);
    
    bDataStyle.dataStyle().fillStyle().setParameter("color","red")
    cDataStyle.dataStyle().fillStyle().setParameter("color","green")
    lDataStyle.dataStyle().fillStyle().setParameter("color","blue")
    bDataStyle.dataStyle().fillStyle().setVisible(true) ;  
    cDataStyle.dataStyle().fillStyle().setVisible(true);
    lDataStyle.dataStyle().fillStyle().setVisible(true);
    
    bDataStyle.dataStyle().errorBarStyle().setVisible(false)
    cDataStyle.dataStyle().errorBarStyle().setVisible(false)
    lDataStyle.dataStyle().errorBarStyle().setVisible(false)
    
    regionStyle = plotter.region(0).style();
    regionStyle.xAxisStyle().setLabel(htitle);
    plotter.setTitle("LCFI Flavour Tag Inputs: " + hname)
    plotter.region(0).setTitle("LCFI Flavour Tag Inputs: " + hname);

    regionStyle.xAxisStyle().labelStyle().setFontSize(16); 
    regionStyle.yAxisStyle().labelStyle().setFontSize(16);
    # can set this next option to "logarithmic"
    regionStyle.yAxisStyle().setScaling("linear");
    regionStyle.titleStyle().textStyle().setFontSize(16);
    regionStyle.yAxisStyle().setLabel("Arbitary Units");

    regionStyle.yAxisStyle().setParameter("allowZeroSuppression","0");
    regionStyle.yAxisStyle().tickLabelStyle().setFontSize(14);  
    regionStyle.xAxisStyle().tickLabelStyle().setFontSize(14); 

    regionStyle.statisticsBoxStyle().setVisibileStatistics("111000");
    regionStyle.statisticsBoxStyle().setVisible(false);
    regionStyle.statisticsBoxStyle().textStyle().setFontSize(14);
   
    regionStyle.legendBoxStyle().textStyle().setFontSize(14); 

    plotter.region(0).plot(bclSum,bDataStyle);
    plotter.region(0).plot(clSum,cDataStyle);
    plotter.region(0).plot(lSum,lDataStyle);

    plotter.show();
    plotter.writeToFile(hname+".png","PNG");   
    # can choose a difference file format - or print both
    #plotter.writeToFile(hname+".eps","EPS");

    return 1;



FlavourTagInputsOverlay("D0Significance1 (zoomed)","Impact parameter significance of most significant track");
FlavourTagInputsOverlay("D0Significance2 (zoomed)","Impact parameter significance of second-most significant track");
#FlavourTagInputsOverlay("DecayLength","Decay Length of most significant track (cm)");
#FlavourTagInputsOverlay("DecayLength(SeedToIP)","Decay Length from Interaction Point (cm)");
#FlavourTagInputsOverlay("JointProbRPhi","Probability of all tracks in jet to come from IP (in r-phi direction)");
#FlavourTagInputsOverlay("JointProbZ","Probability of all tracks in jet to come from IP (in z direction)");
#FlavourTagInputsOverlay("Momentum1","Momentum of most significant track in jet (GeV/c)");
#FlavourTagInputsOverlay("Momentum2","Momentum of second-most significant track in jet (GeV/c)");
#FlavourTagInputsOverlay("NumTracksInVertices","Number of Tracks in non-primary verticies");
#FlavourTagInputsOverlay("NumVertices","Number of vertices per Jet");
#FlavourTagInputsOverlay("PTMassCorrection","PT-Corrected mass of the vertex (GeV/c2)");
#FlavourTagInputsOverlay("RawMomentum","Raw Vertex Momentum (GeV/c)");
#FlavourTagInputsOverlay("SecondaryVertexProbability","Secondary Vertex Probability");
#FlavourTagInputsOverlay("Z0Significance1 (zoomed)","Impact Parameter Significance in z for most significant track in jet");
#FlavourTagInputsOverlay("Z0Significance2 (zoomed)","Impact Parameter Significance in z for second-most significant track in jet");
