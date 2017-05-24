#include <TChain.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <iostream>
#include <RooRealVar.h>
#include <RooVoigtian.h>
#include <RooDataSet.h>
#include <RooChebychev.h>
#include <RooAddPdf.h>
#include <sstream>
#include <TLegend.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <TAxis.h>
#include <TApplication.h>
#include <RooArgList.h>
#include "argparser.h"

using namespace RooFit;
using namespace std;

RooChebychev *CreateChebychev(const char *name, const char *title, const char *parameter_name, RooRealVar &x,
                              vector<double> parameters) {
    RooArgList *p = new RooArgList();
    for (int i = 0; i < parameters.size(); i++) {
        string name(parameter_name);
        name += i;
        string title = "coefficient c#";
        title += i;
        p->add(*(new RooRealVar(name.c_str(), title.c_str(), parameters[i])));
    }
    return new RooChebychev(name, title, x, *p);
}

void Mixkpi_fit_1(string &input_file) {

    TStopwatch timer;
    gSystem->Load("libRooFit");   //load roofit libraries
    TFile *f1 = new TFile(input_file.c_str(), "READ");
    TTree *t1 = (TTree *) f1->Get("treeana");

//####################################for K+ Pi0##################################################
    RooRealVar m_kpstar("m_kpstar", "m_kpstar", 0.62, 1.12);
    RooDataSet datakp("datakp", "Dataset of m_kpstar", t1, m_kpstar);
/////////// Probability Distribution function_1 from K^{+*}K^{-} of f_{1}(1420) decay/////////////
    RooRealVar mass1("mass1", "mass1", 0.8866);
    RooRealVar vsigma1("vsigma1", "vsigma1", 0.026);
    RooRealVar width1("width1", "width1", 0.0195);
    RooVoigtian kppdf1("kppdf1", "from K^{+*}K^{-} of f_{1}(1420) decay", m_kpstar, mass1, vsigma1, width1);
//////////// Probability Distribution function_2 from K^{-*}K^{+} of f_{1}(1420) decay//////////////
//    auto kppdf2 = *(CreateChebychev("kppdf2", "from K^{-*}K^{+} of f_{1}(1420) decay", m_kpstar,
//                                    vector<double>{-0.67, -0.946, 1, -0.082, -0.422, 0.28}));
/////////// Probability Distribution function_3 from K^{+}K^{-}#pi^{0} of f_{1}(1420) decay///////////
//    RooRealVar c_1("c_1", "coefficient c#1", -0.790);
//    RooRealVar c_2("c_2", "coefficient c#2", -0.750);
//    RooRealVar c_3("c_3", "coefficient c#3", 0.868);
//    RooRealVar c_4("c_4", "coefficient c#4", -0.247);
//    RooChebychev kppdf3("kppdf3", "from K^{+}K^{-}#pi^{0} of f_{1}(1420) decay", m_kpstar,
//                        RooArgList(c_1, c_2, c_3, c_4));
/////////// Probability Distribution function_4 from K^{+*}K^{-} of J/psi threebody decay///////////
    RooRealVar mass3("mass3", "mass3", 0.8920);
    RooRealVar vsigma3("vsigma3", "vsigma3", 0.049);
    RooRealVar width3("width3", "width3", 0.008);
    RooVoigtian kppdf4("kppdf4", "from K^{+*}K^{-} of J/psi threebody decay", m_kpstar, mass3, vsigma3, width3);
/////////// Probability Distribution function_5 from K^{-*}K^{+} of J/psi threebody decay///////////
    auto kppdf5 = *(CreateChebychev("kppdf5", "Background", "cc", m_kpstar, vector<double>{-0.28, -1, -0.47, -0.604}));
    /////////// Probability Distribution function_6 from K^{+*}K^{-} of #eta(1475) decay/////////////
    RooRealVar mass5("mass5", "mass5", 0.8779);
    RooRealVar vsigma5("vsigma5", "vsigma5", 0.035);
    RooRealVar width5("width5", "width5", 0.0189);
    RooVoigtian kppdf6("kppdf6", "from K^{+*}K^{-} of #eta(1475) decay", m_kpstar, mass5, vsigma5, width5);
    //////////// Probability Distribution function_7 from K^{-*}K^{+} of #eta(1475) decay//////////////
//    RooRealVar cc1_1("cc1_1", "coefficient c#1", -0.08);
//    RooRealVar cc1_2("cc1_2", "coefficient c#2", -0.276);
//    RooRealVar cc1_3("cc1_3", "coefficient c#3", 0.910);
//    RooRealVar cc1_4("cc1_4", "coefficient c#4", -0.023);
//    RooChebychev kppdf7("kppdf7", "from K^{-*}K^{+} of #eta(1475) decay", m_kpstar,
//                        RooArgList(cc1_1, cc1_2, cc1_3, cc1_4));
//    auto kppdf7 = *(CreateChebychev("kppdf7", "from K^{-*}K^{+} of #eta(1475) decay","cc1",m_kpstar,
//                                    vector<double>{-0.08, -0.276, 0.91, -0.023}));
    /////////// Probability Distribution function_8 from K^{+}K^{-}#pi^{0} of #eta(1475) decay///////////
    auto kppdf8 = *(CreateChebychev("kppdf8", "from K^{+}K^{-}#pi^{0} of #eta(1475) decay", "cc2", m_kpstar,
                                    vector<double>{-0.361, -0.306, 0.87, -0.093}));

///////////////////////////////////sum ////////////////////////////////////////////////////////
    RooRealVar frac1("frac1", " from K^{+*}K^{-} of f_{1}(1420) decay fraction", 0.5, 0, 0.7);
    RooRealVar frac2("frac2", " from K^{-*}K^{+} of f_{1}(1420) decay fraction", 0.3, 0, 0.7);
    RooRealVar frac3("frac3", " from K^{+}K^{-}#pi^{0} of f_{1}(1420) decay fraction", 0.1, 0, 0.8);
    RooRealVar frac4("frac4", " from K^{+*}K^{-} of J/psi threebody decay fraction", 0.1, 0, 0.8);
    RooRealVar frac5("frac5", " from K^{-*}K^{+} of J/psi threebody decay fraction", 0.1, 0, 0.8);
    RooRealVar frac6("frac6", " from K*+ K- of #eta(1475) decay fraction", 0.1, 0, 0.8);
    RooRealVar frac7("frac7", " from K*- K+ of #eta(1475) decay fraction", 0.1, 0, 0.8);
    RooRealVar frac8("frac8", " from K+ K- #pi^{0} of #eta(1475) decay fraction", 0.1, 0, 0.8);
    RooAddPdf kppdf("kppdf", "S1+S2+S3+S4+S5", RooArgList(kppdf4, kppdf5, kppdf8), RooArgList(frac1, frac2));
/////////////////////////////fit//////////////////////////////////////////////////////////////
    auto rf1 = kppdf.fitTo(datakp, Save(true), Strategy(1), NumCPU(8));
    rf1->Print("v");
////////////////////////////////////Plot///////////////////////////////////////////////////////and J/psi 3-body(#omegaK^{+*}K^{-}) decay  #eta(1475) and J/psi 3-body(#omegaK^{#pm*}K^{#mp})
    TString legname1[5];
    RooPlot *xframe1 = m_kpstar.frame(
            Title("the mass of K^{+} #pi^{0} from f_{1}(1420) #eta(1475) and J/psi 3-body(#omegaK^{#pm*}K^{#mp}) decay(K^{+*}K^{-}:K^{-*}K^{+}:K^{+}K^{-}#pi^{0}= 10:10:1)"),
            Bins((1.12 - 0.62) / 0.005));
    datakp.plotOn(xframe1);
    legname1[0] = "data";
    kppdf.plotOn(xframe1, LineColor(1));
    legname1[1] = "fit";
    kppdf.plotOn(xframe1, Components("kppdf4"), LineStyle(kDashed), LineColor(2));
    legname1[2] = "from K^{*+}K^{-} ";
    kppdf.plotOn(xframe1, Components("kppdf5"), LineStyle(kDashed), LineColor(3));
    legname1[3] = "from K^{*-}K^{+} ";
    kppdf.plotOn(xframe1, Components("kppdf8"), LineStyle(kDashed), LineColor(4));
    legname1[4] = "from K^{+}K^{-}#pi^{0} ";


    kppdf.paramOn(xframe1, Format("NEU", AutoPrecision(1)), Layout(0.60, 0.90, 0.90),
                  Parameters(RooArgSet(frac1, frac2)));
//	datakp.statOn(xframe1,Layout(0.1,0.35,0.66));

    ostringstream chi2string1;
    Int_t nParsToFit1 = (rf1->floatParsFinal()).getSize();
    Double_t chi2 = xframe1->chiSquare(nParsToFit1);
    chi2string1 << "#chi^{2}/ndof=" << chi2;
    cout << "#chi^{2}/ndof=" << chi2 << endl;

    TCanvas *c = new TCanvas("c", "k+ pi0 fit plot", 800, 600);
    c->cd();
    xframe1->GetXaxis()->SetTitle("m_{K^{+}#pi^{0}}");
    xframe1->Draw();
    TLegend *legend1 = new TLegend(0.67, 0.52, 0.90, 0.77);
    for (unsigned int i = 0; i < 5; i++) {
        TString objName = xframe1->nameOf(i);
        TObject *obj = xframe1->findObject(objName.Data());
        legend1->AddEntry(obj, legname1[i], "PL");
        legend1->SetTextFont(42);
    }
    legend1->Draw();
//	c->SaveAs("mixkpi_1420_kp_1.pdf");
//####################################for K- Pi0##################################################
    RooRealVar m_kmstar("m_kmstar", "m_kmstar", 0.62, 1.12);
    RooDataSet datakm("datakm", "Dataset of m_kmstar", t1, m_kmstar);
////////////////////////////////// Probability Distribution function_1 /////////////////////////////
    RooRealVar c1_1("c1_1", "coefficient c#1", -0.642);
    RooRealVar c1_2("c1_2", "coefficient c#2", -0.937);
    RooRealVar c1_3("c1_3", "coefficient c#3", 1);
    RooRealVar c1_4("c1_4", "coefficient c#4", -0.057);
    RooRealVar c1_5("c1_5", "coefficient c#5", -0.416);
    RooRealVar c1_6("c1_6", "coefficient c#6", 0.28);
    RooChebychev kmpdf1("kmpdf1", "from K^{+*}K^{-} of f_{1}(1420) decay", m_kmstar,
                        RooArgList(c1_1, c1_2, c1_3, c1_4, c1_5, c1_6));
    ////////////////////////////////// Probability Distribution function_2 /////////////////////////////
    RooRealVar mass2("mass2", "mass2", 0.8870);
    RooRealVar vsigma2("vsigma2", "vsigma2", 0.0251);
    RooRealVar width2("width2", "width2", 0.0191);
    RooVoigtian kmpdf2("kmpdf2", "from K^{-*}K^{+} of f_{1}(1420) decay", m_kmstar, mass2, vsigma2, width2);
    ////////////////////////////////// Probability Distribution function_3 /////////////////////////////
    RooRealVar c_km1("ckm1", "coefficient c#1", -0.796);
    RooRealVar c_km2("ckm2", "coefficient c#2", -0.740);
    RooRealVar c_km3("ckm3", "coefficient c#3", 0.88);
    RooRealVar c_km4("ckm4", "coefficient c#4", -0.260);
    RooChebychev kmpdf3("kmpdf3", "from K^{+}K^{-}#pi^{0} of f_{1}(1420) decay", m_kmstar,
                        RooArgList(c_km1, c_km2, c_km3, c_km4));
    /////////// Probability Distribution function_4 from K^{+*}K^{-} of J/psi threebody decay///////////
    RooRealVar c2_1("c2_1", "coefficient c#1", -0.206);
    RooRealVar c2_2("c2_2", "coefficient c#2", -1);
    RooRealVar c2_3("c2_3", "coefficient c#3", -0.558);
    RooRealVar c2_4("c2_4", "coefficient c#4", -0.653);
    RooChebychev kmpdf4("kmpdf4", "from K^{+*}K^{-} of J/psi threebody decay", m_kmstar,
                        RooArgList(c2_1, c2_2, c2_3, c2_4));
    /////////// Probability Distribution function_5 from K^{-*}K^{+} of J/psi threebody decay///////////
    RooRealVar mass4("mass4", "mass4", 0.8915);
    RooRealVar vsigma4("vsigma4", "vsigma4", 0.049);
    RooRealVar width4("width4", "width4", 0.009);
    RooVoigtian kmpdf5("kmpdf5", "from K^{-*}K^{+} of J/psi threebody decay", m_kmstar, mass4, vsigma4, width4);
    /////////////////// Probability Distribution function_6 from K^{+*}K^{-} of #eta(1475) decay////////////
    RooRealVar c3_1("c3_1", "coefficient c#1", -0.12);
    RooRealVar c3_2("c3_2", "coefficient c#2", -0.392);
    RooRealVar c3_3("c3_3", "coefficient c#3", 0.965);
    RooRealVar c3_4("c3_4", "coefficient c#4", -0.084);
    RooChebychev kmpdf6("kmpdf6", "from K^{+*}K^{-} of #eta(1475) decay", m_kmstar, RooArgList(c3_1, c3_2, c3_3, c3_4));
    /////////////////// Probability Distribution function_7 from K^{-*}K^{+} of #eta(1475) decay////////////
    RooRealVar mass6("mass6", "mass6", 0.8785);
    RooRealVar vsigma6("vsigma6", "vsigma6", 0.035);
    RooRealVar width6("width6", "width6", 0.0185);
    RooVoigtian kmpdf7("kmpdf7", "from K^{-*}K^{+} of #eta(1475) decay", m_kmstar, mass6, vsigma6, width6);
    ////////////////// Probability Distribution function_8 from K^{+}K^{-}#pi^{0} of #eta(1475) decay///////
    RooRealVar c4_1("c4_1", "coefficient c#1", -0.464);
    RooRealVar c4_2("c4_2", "coefficient c#2", -0.412);
    RooRealVar c4_3("c4_3", "coefficient c#3", 0.93);
    RooRealVar c4_4("c4_4", "coefficient c#4", -0.144);
    RooChebychev kmpdf8("kmpdf8", "Background", m_kmstar, RooArgList(c4_1, c4_2, c4_3, c4_4));

///////////////////////////////////sum //////////////////////////////////////////////////////////////
    RooRealVar frc1("frc1", "from K^{+*}K^{-} of f_{1}(1420) decay fraction", 0.5, 0, 0.6);
    RooRealVar frc2("frc2", "from K^{-*}K^{+} of f_{1}(1420) decay fraction", 0.25, 0, 0.6);
    RooRealVar frc3("frc3", "from K^{+}K^{-}#pi^{0} of f_{1}(1420) decay fraction", 0.1, 0, 0.8);
    RooRealVar frc4("frc4", "from K^{+*}K^{-} of J/psi threebody decay fraction", 0.1, 0, 0.8);
    RooRealVar frc5("frc5", "from K^{-*}K^{+} of J/psi threebody decay fraction", 0.2, 0, 0.8);
    RooRealVar frc6("frc6", " from K*+K- of #eta(1475) decay fraction", 0.1, 0, 0.8);
    RooRealVar frc7("frc7", " from K*-K+ of #eta(1475) decay fraction", 0.1, 0, 0.8);
    RooRealVar frc8("frc8", " from K+K-#pi^{0} of #eta(1475) decay fraction", 0.1, 0, 0.8);
    RooAddPdf kmpdf("kmpdf", "S1+S2+S3", RooArgList(kmpdf4, kmpdf5, kmpdf8), RooArgList(frc1, frc2));
/////////////////////////////fit//////////////////////////////////////////////////////////////
    RooFitResult *rf2;
    rf2 = kmpdf.fitTo(datakm, Extended(1), Save(true), Strategy(1), NumCPU(8));
    rf2->Print("v");
////////////////////////////////////Plot///////////////////////////////////////////////////////and J/psi 3-body(#omegaK^{+*}K^{-})  #eta(1475) and J/psi 3-body(#omegaK^{#pm*}K^{#mp})
    TString legname2[5];
    RooPlot *xframe2 = m_kmstar.frame(
            Title("the mass of K^{-} #pi^{0} from f_{1}(1420) #eta(1475) and J/psi 3-body(#omegaK^{#pm*}K^{#mp}) decay(K^{+*}K^{-}:K^{-*}K^{+}:K^{+}K^{-}#pi^{0} = 10:10:1)"),
            Bins((1.12 - 0.62) / 0.005));
    datakm.plotOn(xframe2);
    legname2[0] = "data";
    kmpdf.plotOn(xframe2, LineColor(1));
    legname2[1] = "fit";
    kmpdf.plotOn(xframe2, Components("kmpdf4"), LineStyle(kDashed), LineColor(2));
    legname2[2] = "from K^{*+}K^{-} ";
    kmpdf.plotOn(xframe2, Components("kmpdf5"), LineStyle(kDashed), LineColor(3));
    legname2[3] = "from K^{*-}K^{+} ";
    kmpdf.plotOn(xframe2, Components("kmpdf8"), LineStyle(kDashed), LineColor(4));
    legname2[4] = "from K^{+}K^{-}#pi^{0} ";


    kmpdf.paramOn(xframe2, Format("NEU", AutoPrecision(1)), Layout(0.60, 0.90, 0.90),
                  Parameters(RooArgSet(frc1, frc2)));
//	datakm.statOn(xframe2,Layout(0.1,0.35,0.66));

    ostringstream chi2string2;
    Int_t nParsToFit2 = (rf2->floatParsFinal()).getSize();
    chi2 = xframe2->chiSquare(nParsToFit2);
    chi2string2 << "#chisq^{2}/ndof=" << chi2;
    cout << "#chisq^{2}/ndof=" << chi2 << endl;

    TCanvas *cc = new TCanvas("cc", "k-pi0 fit plot", 800, 600);
    cc->cd();
    xframe2->GetXaxis()->SetTitle("m_{K^{-}#pi^{0}}");
    xframe2->Draw();
    TLegend *legend2 = new TLegend(0.67, 0.52, 0.90, 0.77);
    for (unsigned int i = 0; i < 5; i++) {
        TString objName = xframe2->nameOf(i);
        TObject *obj = xframe2->findObject(objName.Data());
        legend2->AddEntry(obj, legname2[i], "PL");
        legend2->SetTextFont(42);
    }
    legend2->Draw();

//	cc->SaveAs("mixkpi_1420_km_1.pdf");
    ///////////////////////////////estimate fit time////////////////////////////////////
    double cputime = timer.CpuTime();
    printf("RT=%7.3f s, Cpu=%7.3f s\n", timer.RealTime(), cputime);

}

int main(int argc, char *argv[]) {
    ArgumentParser parser;
    parser.addArgument("-i", "--input", 1, false);
    parser.parse(argc, argv);
    TApplication *myapp = new TApplication("App", &argc, argv);
    Mixkpi_fit_1(parser.retrieve<string>("input"));
    myapp->Run();
}