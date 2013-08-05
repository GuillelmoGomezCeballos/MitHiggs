void MakeMitClass()
{
  TFile* infile = new TFile("/home/ceballos/condor/old/histo_f10-wz-z2-v12_all_noskim.root","read");

  assert(infile);

  HwwTree->MakeClass("MitNtupleEvent");
}
