//--------------------------------------------------------------------------------------------------
// $Id $
//
// PlotUtils
//
// Utility functions for plotting
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITHIGGS_UTILS_PLOTUTILS_H
#define MITHIGGS_UTILS_PLOTUTILS_H

#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/Utils/interface/SimpleTable.h"
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include "MitCommon/DataFormats/interface/TH3DAsymErr.h"
#include <string>
#include <vector>
#include <utility>

namespace mithep
{
  class PlotUtils {

    public:
      PlotUtils();
      PlotUtils(Double_t x);
      ~PlotUtils();

      void Init(Double_t x);
      void SetLuminosity(Double_t x)         { fIntegratedLuminosity = x ; }
      std::vector<Int_t> GetCOLORS()         { return fCOLORS;             }
      std::vector<Int_t> GetSYSCOLORS()      { return fSYSCOLORS;          }
      std::vector<Int_t> GetMARKERS()        { return fMARKERS;            }

      Double_t getWeight(std::string datasetFile, std::string datasetName);
      static TH1F *getHisto(std::string filename, std::string directoryname, std::string histoname);
      static TH2F* get2DHisto(std::string filename, std::string directoryname, std::string histoname);
      static TH1F* rebin(TH1F* hist, std::vector<Double_t> xlowedges, const string name = "");
      static TH1F* rebin(TH1F* hist, Int_t nbins, const string name = "");
      static TH2F* rebin(TH2F* hist, std::vector<Double_t> xlowedges, std::vector<Double_t> ylowedges, const string name = "");
      static TH3F* rebin(TH3F* hist, std::vector<Double_t> xlowedges, std::vector<Double_t> ylowedges, 
                  std::vector<Double_t> zlowedges);
      TH1F* addAllSamples(std::vector<std::string> datasetFiles, 
                          std::vector<std::string> datasetNames,
                          std::string dirName, std::string histName, 
                          std::vector<Double_t> bins, 
                          std::string xAxisTitle, std::string yAxisTitle);
      TH2F* addAllSamples2D(std::vector<std::string> datasetFiles, 
                            std::vector<std::string> datasetNames,
                            std::string dirName, std::string histName, 
                            std::vector<Double_t> xbins, 
                            std::vector<Double_t> ybins);
      THStack* addAllSamplesStacked(std::vector<std::string> datasetFiles, 
                                    std::vector<std::string> datasetNames, 
                                    std::string dirName, std::string histName, 
                                    std::vector<Double_t> bins, 
                                    std::string xAxisTitle,
                                    std::string yAxisTitle);
      THStack* addAllSamplesStacked(std::vector<std::vector<std::string> > datasetFiles, 
                                    std::vector<std::vector<std::string> > datasetNames, 
                                    std::string dirName, std::string histName, 
                                    std::vector<Double_t> bins, 
                                    std::string xAxisTitle,
                                    std::string yAxisTitle);
      TH1F* createEfficiencyHist(std::vector<std::string> datasetFiles, 
                                 std::vector<std::string> datasetNames, 
                                 std::string dirName,std::string numeratorHistname, 
                                 std::string denominatorHistname, 
                                 std::string histname, std::vector<Double_t> bins, 
                                 Double_t xlow, Double_t xhigh, 
                                 Double_t ylow, Double_t yhigh);
      TGraphAsymmErrors* createEfficiencyGraph(std::vector<std::string> datasetFiles, 
                                             std::vector<std::string> datasetNames, 
                                             std::string dirName,std::string numeratorHistname, 
                                             std::string denominatorHistname, 
                                             std::string histname, std::vector<Double_t> bins, 
                                             Int_t errorType,
                                             Double_t xlow, Double_t xhigh, 
                                             Double_t ylow, Double_t yhigh);
      TGraphAsymmErrors* createEfficiencyGraph(TH1F* numerator, TH1F* denominator,
                                               std::string histname, std::vector<Double_t> bins, 
                                               Int_t errorType,
                                               Double_t xlow, Double_t xhigh, 
                                               Double_t ylow, Double_t yhigh);
     

      void drawStackedPlot(THStack *stackedHist , std::string plotname, 
                                  std::vector<std::string> legendNames,
                                  Bool_t logY, Double_t MaxY, Double_t MinX, Double_t MaxX,
                                  Double_t legendX1, Double_t legendY1, 
                                  Double_t legendX2, Double_t legendY2,
                                  TH1F *hist, std::string histLegendLabel);
      TCanvas* DrawDataSignalBkgHistogram(TH1F* data, TH1F* sig, THStack *bkg , TLegend *legend, 
                                          std::string plotname, 
                                          Bool_t useLogY, double MaxY, 
                                          double MinX, double MaxX,
                                          double legendX1, double legendY1, 
                                          double legendX2, double legendY2);
      void makeDistributionComparisonPlot( std::vector<std::string> datasetfiles, 
                                           std::vector<std::string> datasetnames, 
                                           std::string dirName,
                                           std::vector<std::string> histNames, 
                                           std::vector<std::string> legendNames, 
                                           bool normalizeArea,
                                           Double_t MinX, Double_t MaxY,
                                           int nbins, std::string plotname );
      void makeXSliceDistributionComparisonPlot( std::vector<TH2F*> hists,                                       
                                                 std::vector<std::string> legendNames, 
                                                 bool normalizeArea,
                                                 std::string xAxisLabel,
                                                 std::string yAxisLabel,
                                                 Double_t xlow, Double_t xhigh, 
                                                 Double_t ylow, Double_t yhigh, 
                                                 double legendX1, double legendX2 , 
                                                 double legendY1, double legendY2,
                                                 std::string plotname );
        void makeXSliceDistributionComparisonPlot( std::vector<std::string> datasetfiles, 
                                                   std::vector<std::string> datasetnames, 
                                                   std::string dirName,
                                                   std::vector<std::string> histNames, 
                                                   std::vector<std::string> legendNames, 
                                                   bool normalizeArea,
                                                   std::vector<double> xbins, std::vector<double> ybins,
                                                   std::string xAxisLabel,
                                                   std::string yAxisLabel,
                                                   Double_t xlow, Double_t xhigh, 
                                                   Double_t ylow, Double_t yhigh, 
                                                   double legendX1, double legendX2 , 
                                                   double legendY1, double legendY2, 
                                                   std::string plotname );
      void makeYSliceDistributionComparisonPlot( std::vector<TH2F*> hists,     
                                                 std::vector<std::string> legendNames, 
                                                 bool normalizeArea,
                                                 std::string xAxisLabel,
                                                 std::string yAxisLabel,
                                                 Double_t xlow, Double_t xhigh, 
                                                 Double_t ylow, Double_t yhigh, 
                                                 double legendX1, double legendX2 , 
                                                 double legendY1, double legendY2,
                                                 std::string plotname );
      void makeYSliceDistributionComparisonPlot( std::vector<std::string> datasetfiles, 
                                                 std::vector<std::string> datasetnames, 
                                                 std::string dirName,
                                                 std::vector<std::string> histNames, 
                                                 std::vector<std::string> legendNames, 
                                                 bool normalizeArea,
                                                 std::vector<double> xbins, std::vector<double> ybins,
                                                 std::string xAxisLabel,
                                                 std::string yAxisLabel,
                                                 Double_t xlow, Double_t xhigh, 
                                                 Double_t ylow, Double_t yhigh, 
                                                 double legendX1, double legendX2 , 
                                                 double legendY1, double legendY2,
                                                 std::string plotname );
      void makeCrossDatasetComparisonPlot( std::vector<std::string> dataset1files, 
                                           std::vector<std::string> dataset1names, 
                                           std::string dataset1label,
                                           std::vector<std::string> dataset2files, 
                                           std::vector<std::string> dataset2names, 
                                           std::string dataset2label,
                                           std::string dirname1,
                                           std::vector<std::string> histNames1, 
                                           std::vector<std::string> legendNames1,
                                           std::string dirname2,
                                           std::vector<std::string> histNames2, 
                                           std::vector<std::string> legendNames2, 
                                           std::string plotname,
                                           bool normalizeArea, 
                                           std::string xAxisLabel,
                                           std::string yAxisLabel,
                                           Double_t xlow, Double_t xhigh, 
                                           Double_t ylow, Double_t yhigh, 
                                           Double_t legendX1, Double_t legendX2, 
                                           Double_t legendY1, Double_t legendY2, 
                                           Bool_t useLogY, Int_t nbins);
      void makeCrossDirComparisonPlot( std::vector<std::string> datasetfiles, 
                                       std::vector<std::string> datasetnames, 
                                       std::vector<std::string> dirnames, 
                                       std::vector<std::string> dirnamelabel, 
                                       std::vector<std::string> histNames,
                                       std::vector<std::string> legendNames, 
                                       std::vector<Double_t> bins, std::string plotname , 
                                       std::string xAxisLabel,
                                       std::string yAxisLabel,
                                       Double_t xlow, Double_t xhigh, 
                                       Double_t ylow, Double_t yhigh, 
                                       Double_t legendX1, Double_t legendX2, 
                                       Double_t legendY1, Double_t legendY2, 
                                       Bool_t useLogY );
      void makeComparisonPlotWithSystematics( std::vector<std::string> datasetfiles, 
                                              std::vector<std::string> datasetnames, 
                                              std::vector<std::string> dirnames, 
                                              std::vector<std::string> dirnamelabel, 
                                              std::vector<std::string> histNames, 
                                              std::vector<std::string> errorhistNames, 
                                              std::vector<std::string> legendNames, 
                                              std::vector<Double_t> bins, 
                                              std::string plotname,
                                              std::string xAxisLabel,
                                              std::string yAxisLabel,
                                              Double_t xlow, Double_t xhigh, 
                                              Double_t ylow, Double_t yhigh, 
                                              Double_t legendX1, Double_t legendX2, 
                                              Double_t legendY1, Double_t legendY2, 
                                              Bool_t useLogY);
   
      static void NormalizeHist(TH1F *hist);
      static TGraphAsymmErrors* MakeSigEffVsBkgEffGraph(TH1F* signalHist, TH1F* bkgHist, 
                                                        std::string name, 
                                                        Bool_t cutBelow = kTRUE );
      static TGraphAsymmErrors* MakeCurrentWPSigEffVsBkgEffGraph(Double_t signalEff, Double_t bkgEff, 
                                                                 std::string name );
      static TGraphAsymmErrors* MakeCurrentWPSigEffVsBkgEffGraph(TH1F* signalHist, TH1F* bkgHist, 
                                                                 std::string name, Double_t myCutValue, 
                                                                 Bool_t cutBelow = kTRUE);
      static TGraphAsymmErrors* MakeSigEffVsCutValueGraph(TH1F* signalHist, std::string name , 
                                                          Bool_t cutBelow = kTRUE);
      static TGraphAsymmErrors* MakeCurrentWPSigEffVsCutValueGraph(TH1F* signalHist, std::string name, 
                                                                   Double_t myCutValue, 
                                                                   Bool_t cutBelow  = kTRUE);
      static Double_t FindCutValueAtFixedSignalEfficiency(TH1F* signalHist, Double_t targetSignalEff, 
                                                          Bool_t cutBelow = kTRUE );
      static Double_t FindBkgEffAtFixedSignalEfficiency(TH1F* signalHist, TH1F* bkgHist, 
                                                        Double_t targetSignalEff, 
                                                        Bool_t cutBelow = kTRUE );
      

    private:
      std::vector<Int_t> fCOLORS;
      std::vector<Int_t> fMARKERS; 
      std::vector<Int_t> fSYSCOLORS;
      Double_t fIntegratedLuminosity;
      SimpleTable *xstab;

      ClassDef(PlotUtils, 0) // Plot utility functions
  };
}

#endif
