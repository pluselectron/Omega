#include <TSystem.h>
#include <TStopwatch.h>
#include <TFile.h>
#include <TTree.h>
#include <RooRealVar.h>
#include <RooKeysPdf.h>
#include <RooDataSet.h>
//#include <TApplication.h>

//
// Created by 郑逸 on 17/5/24.
//
void getpdf() {
    TStopwatch timer;
    gSystem->Load("libRooFit");   //load roofit libraries

    TFile *f6 = new TFile("~/workshop/omegaG2K2/Ana_omegaG2K2_2/signalMC/test_ShapeFit_method/1475/omegaEta1475_kp_1000wan_cutWithoutKstar.root");
    TTree *t6 = (TTree *) f6->Get("treeana");
    TFile *f7 = new TFile("~/workshop/omegaG2K2/Ana_omegaG2K2_2/signalMC/test_ShapeFit_method/1475/omegaEta1475_km_1000wan_cutWithoutKstar.root");
    TTree *t7 = (TTree *) f7->Get("treeana");
    TFile *f8 = new TFile("~/workshop/omegaG2K2/Ana_omegaG2K2_2/signalMC/test_ShapeFit_method/1475/omegaEta1475_dir_100wan_cutWithoutKstar.root");
    TTree *t8 = (TTree *) f8->Get("treeana");
    /*
    TFile *f6 =new TFile("omegaEta1475_kps.root");
    TTree *t6 = (TTree*)f6->Get("treeana");
    TFile *f7 = new TFile("omegaEta1475_kms.root");
    TTree *t7 = (TTree*)f7->Get("treeana");
    TFile *f8 = new TFile("omegaEta1475_dir.root");
    TTree *t8 = (TTree*)f8->Get("treeana");
  */
    //The range of x
    RooRealVar m_kpstar("m_kpstar", "m_kpstar", 0.62, 1.62);
    RooRealVar m_kmstar("m_kmstar", "m_kmstar", 0.62, 1.62);

    RooDataSet datakp6("datakp6", "K^{+}#pi^{0} from K*+K- of #eta(1475) decay", t6, m_kpstar);
    RooDataSet datakp7("datakp7", "K^{+}#pi^{0} from K*-K+ of #eta(1475) decay", t7, m_kpstar);
    RooDataSet datakp8("datakp8", "K^{+}#pi^{0} from K+K-#pi^{0} of #eta(1475) decay", t8, m_kpstar);

    RooDataSet datakm6("datakm6", "K^{-}#pi^{0} from K*+K- of #eta(1475) decay", t6, m_kmstar);
    RooDataSet datakm7("datakm7", "K^{-}#pi^{0} from K*-K+ of #eta(1475) decay", t7, m_kmstar);
    RooDataSet datakm8("datakm8", "K^{-}#pi^{0} from K+K-#pi^{0} of #eta(1475) decay", t8, m_kmstar);

    RooKeysPdf kppdf6("kppdf6", "K^{+}#pi^{0} from K*+K- of #eta(1475) decay", m_kpstar, datakp6, RooKeysPdf::NoMirror,
                      2);
    RooKeysPdf kppdf7("kppdf7", "K^{+}#pi^{0} from K*-K+ of #eta(1475) decay", m_kpstar, datakp7, RooKeysPdf::NoMirror,
                      2);
    RooKeysPdf kppdf8("kppdf8", "K^{+}#pi^{0} from K+K-#pi^{0} of #eta(1475) decay", m_kpstar, datakp8,
                      RooKeysPdf::NoMirror, 2);

    RooKeysPdf kmpdf6("kmpdf6", "K^{-}#pi^{0} from K*+K- of #eta(1475) decay", m_kmstar, datakm6, RooKeysPdf::NoMirror,
                      2);
    RooKeysPdf kmpdf7("kmpdf7", "K^{-}#pi^{0} from K*-K+ of #eta(1475) decay", m_kmstar, datakm7, RooKeysPdf::NoMirror,
                      2);
    RooKeysPdf kmpdf8("kmpdf8", "K^{-}#pi^{0} from K+K-#pi^{0} of #eta(1475) decay", m_kmstar, datakm8,
                      RooKeysPdf::NoMirror, 2);

    TFile *f = new TFile("./1475/pdf.root", "recreate");
    kppdf6.Write();
    kppdf7.Write();
    kppdf8.Write();
    kmpdf6.Write();
    kmpdf7.Write();
    kmpdf8.Write();
    f->Close();
}

int main(int argc, char *argv[])
{
//    TApplication *myapp = new TApplication("App", &argc, argv);
    getpdf();
//    myapp->Run();
}
