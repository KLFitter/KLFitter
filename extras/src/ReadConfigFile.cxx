#include "ReadConfigFile.h" 
#include <iostream> 
#include <cmath>
#include <fstream>
#include <TSystem.h>  

// --------------------------------------------------------- 
KLFitter::ReadConfigFile::ReadConfigFile(std::string filename)
{
  LeptonType = KLFitter::LikelihoodTopLeptonJets::kElectron;

  BTaggingMethod     = KLFitter::LikelihoodBase::kNotag;
  CutBTagging      = 1.e4;
  FlagIntegrate    = false;
  FlagTopMassFixed = false;
  FlagHiggsMassFixed = false;
  FlagWriteSignalMCTruth   = true;
  BeamEnergy = KLFitter::DetectorBase::k7TeV;
  FlagTruthSel = false;
  LJetSeparationMethod = KLFitter::LikelihoodTopLeptonJetsUDSep::kNone;
  TopPoleMass = 172.5;
  HiggsMass = 120.0;
  
  input_path="input.root";
  output_path="output.root";

  //Following variables are obsolete and not used anymore
  IsBkg = false;
  FlagAthenaComp   = false;
  DO_BATCH    = false;
  FlagUseJetMass   = false;


  ReadConfig(filename);
}


// --------------------------------------------------------- 
KLFitter::ReadConfigFile::ReadConfigFile(std::string filename, bool * validconfig)
{
  LeptonType = KLFitter::LikelihoodTopLeptonJets::kElectron;

  BTaggingMethod     = KLFitter::LikelihoodBase::kNotag;
  CutBTagging      = 1.e4;
  FlagIntegrate    = false;
  FlagTopMassFixed = false;
  FlagHiggsMassFixed = false;
  FlagWriteSignalMCTruth   = true;
  BeamEnergy = KLFitter::DetectorBase::k7TeV;
  FlagTruthSel = false;
  LJetSeparationMethod = KLFitter::LikelihoodTopLeptonJetsUDSep::kNone;

  TopPoleMass = 172.5;
  HiggsMass = 120.0;

  input_path="input.root";
  output_path="output.root";

  //Following variables are obsolete and not used anymore
  IsBkg = false;
  FlagAthenaComp   = false;
  DO_BATCH    = false;
  FlagUseJetMass   = false;

  if(ReadConfig(filename)!=-1){*validconfig=true;}
  else{*validconfig=false;}
}


// --------------------------------------------------------- 
KLFitter::ReadConfigFile::~ReadConfigFile()
{
}

// --------------------------------------------------------- 
int KLFitter::ReadConfigFile::ReadConfig(std::string filename)
{
  ifstream configfile(filename.c_str(),std::ios::in);
  std::string line;
  bool is_comment=false;
  bool is_whitespace=false;
  bool slash=false;
  int tmp;
  double tmpdouble=0.0;
  std::string tmpstr="";
  size_t found;
  std::string::iterator k;
  std::cout<<"reading "<<filename.c_str()<<"..."<<std::endl;

  if(configfile.is_open())
    {
      //read config
      while (! configfile.eof() )
        {
          getline (configfile,line);
          is_comment=false;
          slash=false;
          is_whitespace=false;
          if(line.size()==0) //ignore any empty lines
            continue;
          //looking for uncommentation
          for(k=line.begin();k<line.end();k++)
            {
              //ignore any white space, unless a slash was found before
              if((*k==std::string(" ") || *k==std::string("\t") || *k==std::string("\n")) && slash==false)
                {
                  is_whitespace=true;
                  continue;
                }
              //a line is a comment, if there is a #, % or // as first non whitespace character
              if(*k==std::string("#") || *k==std::string("%") || (*k==std::string("/") && slash==true))
                {
                  is_comment=true;
                  is_whitespace=false;
                  break;
                }
              else if(*k==std::string("/"))
                {
                  slash=true;
                  is_whitespace=false;
                }
              else
                {
                  break;
                }
            }

          if(is_comment==false)
            {
              found=line.find("LeptonType");
              if(found!=std::string::npos)
                {
                  found=line.find("=",found);
                  if(found!=std::string::npos)
                    {
                      tmp=GetPath(&tmpstr,line,found);
                      if(tmp!=-1)
                  {
                    if(tmpstr == "electron")
                      LeptonType=KLFitter::LikelihoodTopLeptonJets::kElectron;
                    if(tmpstr == "muon")
                      LeptonType=KLFitter::LikelihoodTopLeptonJets::kMuon;
                  }
                      else
                        {
                          std::cout<<"Warning: Error while reading value of LeptonType, using standard value"<<std::endl;
                        }
                    }
                }
		  else
		    {
		      found=line.find("CutBTagging");
		      if(found!=std::string::npos)
			{
			  found=line.find("=",found);
			  if(found!=std::string::npos)
			    {
			      tmp=GetValue(&tmpdouble,line,found);
				if(tmp!=-1)
				  {
				    CutBTagging=tmpdouble;
				  }
				else
				  {
				    std::cout<<"Warning: Error while reading value of CutBTagging, using standard value"<<std::endl;
				  }
			    }
			}
		      else
			{
			  found=line.find("DO_BATCH");
			  if(found!=std::string::npos)
			    {
			      found=line.find("=",found);
			      if(found!=std::string::npos)
				{
				  tmp=GetTrueOrFalse(line,found);
				  if(tmp!=-1)
				    {
				      DO_BATCH=(tmp==1);
				    }
				  else
				    {
				      std::cout<<"Warning: Error while reading value of DO_BATCH, using standard value"<<std::endl;
				    }
				}
			    }
			  else
			    {
			      found=line.find("BTaggingMethod");
			      if(found!=std::string::npos)
				{
				  found=line.find("=",found);
          if(found!=std::string::npos)
          {
          tmp=GetPath(&tmpstr,line,found);
          if(tmp!=-1)
          {
          if(tmpstr == "Notag")
            BTaggingMethod=KLFitter::LikelihoodBase::kNotag;
          if(tmpstr == "Veto")
            BTaggingMethod=KLFitter::LikelihoodBase::kVeto;
          if(tmpstr == "VetoLight")
            BTaggingMethod=KLFitter::LikelihoodBase::kVetoLight;
          if(tmpstr == "VetoBoth")
            BTaggingMethod=KLFitter::LikelihoodBase::kVetoBoth;
          if(tmpstr == "VetoNoFit")
            BTaggingMethod=KLFitter::LikelihoodBase::kVetoNoFit;		
          if(tmpstr == "VetoNoFitLight")
            BTaggingMethod=KLFitter::LikelihoodBase::kVetoNoFitLight;	
          if(tmpstr == "VetoNoFitBoth")
            BTaggingMethod=KLFitter::LikelihoodBase::kVetoNoFitBoth;	  
          if(tmpstr == "WorkingPoint")
            BTaggingMethod=KLFitter::LikelihoodBase::kWorkingPoint;
                  }

				      else
					{
					  std::cout<<"Warning: Error while reading value of BTaggingMethod, using standard value"<<std::endl;
					}
				    }
				}
			      else
				{
				  found=line.find("FlagTopMassFixed");
				  if(found!=std::string::npos)
				    {
				      found=line.find("=",found);
				      if(found!=std::string::npos)
					{
					  tmp=GetTrueOrFalse(line,found);
					  if(tmp!=-1){
					    FlagTopMassFixed=(tmp==1);
					  }
					  else
					    {
					      std::cout<<"Warning: Error while reading value of FlagTopMassFixed, using standard value"<<std::endl;
					    }
					}
				    }
				  else
				    {
				      found=line.find("FlagHiggsMassFixed");
				      if(found!=std::string::npos)
					{
					  found=line.find("=",found);
					  if(found!=std::string::npos)
					    {
					      tmp=GetTrueOrFalse(line,found);
					      if(tmp!=-1){
						FlagHiggsMassFixed=(tmp==1);
					      }
					      else
						{
						  std::cout<<"Warning: Error while reading value of FlagHiggsMassFixed, using standard value"<<std::endl;
						}
					    }
					}
				      else
				    {
				      found=line.find("FlagIntegrate");
				      if(found!=std::string::npos)
					{
					  found=line.find("=",found);
					  if(found!=std::string::npos)
					    {
					      tmp=GetTrueOrFalse(line,found);
					      if(tmp!=-1){
						FlagIntegrate=(tmp==1);
					      }
					      else
						{
						  std::cout<<"Warning: Error while reading value of FlagIntegrate, using standard value"<<std::endl;
						}
					    }
					}
				      else
					{
					  found=line.find("FlagUseJetMass");
					  if(found!=std::string::npos)
					    {
					      found=line.find("=",found);
					      if(found!=std::string::npos)
						{
						  tmp=GetTrueOrFalse(line,found);
						  if(tmp!=-1)
						    {
						      FlagUseJetMass=(tmp==1);
						    }
						  else
						    {
						      std::cout<<"Warning: Error while reading value of FlagUseJetMass, using standard value"<<std::endl;
						    }
						}
					    }
					  else
					    {
					      found=line.find("FlagWriteSignalMCTruth");
					      if(found!=std::string::npos)
						{
						  found=line.find("=",found);
						  if(found!=std::string::npos)
						    {
						      tmp=GetTrueOrFalse(line,found);
						      if(tmp!=-1)
							{
							  FlagWriteSignalMCTruth=(tmp==1);
							}
						      else
							{
							  std::cout<<"Warning: Error while reading value of FlagWriteSignalMCTruth, using standard value"<<std::endl;
							}
						    }
						}
					      else
						{
						  found=line.find("BeamCMEnergy");
              if(found!=std::string::npos)
                  {
                    found=line.find("=",found);
                    if(found!=std::string::npos)
                {
                  tmp=GetValue(&tmpdouble,line,found);
                  if(tmp!=-1)
                    {
                      if(tmpdouble==7)
                        BeamEnergy=KLFitter::DetectorBase::k7TeV;
                      if(tmpdouble==8)
                        BeamEnergy=KLFitter::DetectorBase::k8TeV;
                      if(tmpdouble==10)
                        BeamEnergy=KLFitter::DetectorBase::k10TeV;
                    }
                  else
                    {
                      std::cout<<"Warning: Error while reading value of BeamCMEnergy, using standard value"<<std::endl;
                    }
                }
                  }
						  else
							{
							  found=line.find("TopPoleMass");
							  if(found!=std::string::npos)
							    {
							      found=line.find("=",found);
							      if(found!=std::string::npos)
								{
								  tmp=GetValue(&tmpdouble,line,found);
								  if(tmp!=-1)
								    {
								      TopPoleMass=tmpdouble;
								    }
								  else
								    {
								      std::cout<<"Warning: Error while reading value of TopPoleMass, using standard value"<<std::endl;
								    }
								}
							    }
							  else
							    {
							      found=line.find("HiggsMass");
							      if(found!=std::string::npos)
								{
								  found=line.find("=",found);
								  if(found!=std::string::npos)
								    {
								      tmp=GetValue(&tmpdouble,line,found);
								      if(tmp!=-1)
									{
									  HiggsMass=tmpdouble;
									}
								      else
									{
									  std::cout<<"Warning: Error while reading value of HiggsMass, using standard value"<<std::endl;
									}
								    }
								}
							      else
								{
								  found=line.find("PathToInputFile");
								  if(found!=std::string::npos)
								{
								  found=line.find("=",found);
								  if(found!=std::string::npos)
								    {
								      tmp=GetPath(&tmpstr,line,found);
								      if(tmp!=-1)
									{
									  input_path=tmpstr;
									}
								      else
									{
									  std::cout<<"Warning: Error while reading path to input file, using standard path 'input.root'"<<std::endl;
									}
								    }
								}   
							      else
								{
								  found=line.find("PathToOutputFile");
								  if(found!=std::string::npos)
								    {
								      found=line.find("=",found);
								      if(found!=std::string::npos)
									{
									  tmp=GetPath(&tmpstr,line,found);
									  if(tmp!=-1)
									    {
									      output_path=tmpstr;
									    }
									  else
									    {
									      std::cout<<"Warning: Error while reading path to output file, using standard path 'output.root'"<<std::endl;
									    }
									}
								    }
								  else
								    {
								      found=line.find("IsBkg");
								      if(found!=std::string::npos)
									{
									  found=line.find("=",found);
									  if(found!=std::string::npos)
									    {
									      tmp=GetTrueOrFalse(line,found);
									      if(tmp!=-1)
										{
										  IsBkg=(tmp==1);
										}
									      else
										{
										  std::cout<<"Warning: Error while reading value of IsBkg, using standard value"<<std::endl;
										}
									    }
									}
                                                                      else
                                                                        {
                                                                          found=line.find("FlagTruthSel");
                                                                          if(found!=std::string::npos)
                                                                            {
                                                                              found=line.find("=",found);
                                                                              if(found!=std::string::npos)
                                                                                {
                                                                                  tmp=GetTrueOrFalse(line,found);
                                                                                  if(tmp!=-1)
                                                                                    {
                                                                                      FlagTruthSel=(tmp==1);
                                                                                    }
                                                                                  else
                                                                                    {
                                                                                      std::cout<<"Warning: Error while reading value of FlagTruthSel, using standard value"<<std::endl;
                                                                                    }
                                                                                }
                                                                            }
                                                                          else
                                                                            {
                                                                              found=line.find("LJetSeparationMethod");
                                                                              if(found!=std::string::npos)
                                                                                {
                                                                                  found=line.find("=",found);
                                                                                  if(found!=std::string::npos)
                                                                                    {
                                                                                      tmp=GetPath(&tmpstr,line,found);
                                                                                      if(tmp!=-1)
                                                                                        {
                                                                                          if(tmpstr == "None")
											  LJetSeparationMethod=KLFitter::LikelihoodTopLeptonJetsUDSep::kNone;
										          if(tmpstr == "PermReweight")
											  LJetSeparationMethod=KLFitter::LikelihoodTopLeptonJetsUDSep::kPermReweight;
                                                                                        }
                                                                                      else
                                                                                        {
                                                                                          std::cout<<"Warning: Error while reading value of LJetSeparationMethod, using standard value"<<std::endl;
                                                                                        }
                                                                                    }
                                                                                }  
                                                                              else
                                                                              {
                                                                                 found=line.find("FlagAthenaComp");
                                                                                 if(found!=std::string::npos)
                                                                                   {
                                                                                     found=line.find("=",found);
                                                                                     if(found!=std::string::npos)
                                                                                       {
                                                                                         tmp=GetTrueOrFalse(line,found);
                                                                                         if(tmp!=-1)
                                                                                           {
                                                                                             FlagAthenaComp=(tmp==1);
                                                                                           }
                                                                                         else
                                                                                           {
                                                                                             std::cout<<"Warning: Error while reading value of FlagAthenaComp, using standard value"<<std::endl;
                                                                                           }
                                                                                       }
                                                                                   }
                                                                                 else
                                                                                   {
                                                                                     if(is_whitespace==false)
                                                                                       {
                                                                                         std::cout<<"Warning: the line \""<<line.c_str()<<"\" does not match any variable. It is ignored."<<std::endl;
                                                                                       }
                                                                                   }
										}
                                                                            }
									}
								    }
								}
							    }
							}
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
      configfile.close();
    }
  else
    {
      std::cout<<"Couldn't open the file \"";
      std::cout<<filename.c_str();
      std::cout<<"\". Will use standard values."<<std::endl;
    }
  //Print Flag values
  std::cout<<"Flags are set as follows:"<<std::endl;
  if (LeptonType==KLFitter::LikelihoodTopLeptonJets::kElectron)
    std::cout<<"LeptonType = Electron"<<std::endl;
  if (LeptonType==KLFitter::LikelihoodTopLeptonJets::kMuon)
    std::cout<<"LeptonType = Muon"<<std::endl;
  if(BTaggingMethod==KLFitter::LikelihoodBase::kNotag)
    std::cout<< "BTaggingMethod = Notag"<<std::endl;
  if(BTaggingMethod==KLFitter::LikelihoodBase::kVeto)
    std::cout<< "BTaggingMethod = Veto"<<std::endl;
  if(BTaggingMethod==KLFitter::LikelihoodBase::kVetoNoFit)
    std::cout<< "BTaggingMethod = VetoNoFit"<<std::endl;  
  if(BTaggingMethod==KLFitter::LikelihoodBase::kWorkingPoint)
    std::cout<< "BTaggingMethod = WorkingPoint"<<std::endl;
  std::cout<< "CutBTagging = "<<CutBTagging<<std::endl;
  std::cout<< "FlagIntegrate = "<<FlagIntegrate<<std::endl;
  std::cout<< "FlagTopMassFixed = "<<FlagTopMassFixed<<std::endl;
  std::cout<< "FlagHiggsMassFixed = "<<FlagHiggsMassFixed<<std::endl;
  std::cout<< "TopPoleMass = "<<TopPoleMass<<std::endl;
  std::cout<< "HiggsMass = "<<HiggsMass<<std::endl;
  std::cout<< "FlagUseJetMass = "<<FlagUseJetMass<<std::endl;
  std::cout<< "FlagWriteSignalMCTruth = "<<FlagWriteSignalMCTruth<<std::endl;
  std::cout << "FlagTruthSel = " << FlagTruthSel << std::endl;
  if(LJetSeparationMethod==KLFitter::LikelihoodTopLeptonJetsUDSep::kNone)
    std::cout << "LJetSeparationMethod = None" << std::endl;
  if(LJetSeparationMethod==KLFitter::LikelihoodTopLeptonJetsUDSep::kPermReweight)
    std::cout << "LJetSeparationMethod = PermReweight" << std::endl;
  if (BeamEnergy==KLFitter::DetectorBase::k7TeV)
    std::cout<< "BeamCMEnergy = 7TeV" <<std::endl;
  if (BeamEnergy==KLFitter::DetectorBase::k10TeV)
    std::cout<< "BeamCMEnergy = 10TeV" <<std::endl;
  std::cout<< "PathToInputFile = " << input_path<<std::endl;
  std::cout<< "PathToOutputFile = " << output_path<<std::endl;
  std::cout << std::endl;

  if(CheckIOPath()){return 0;}
  else{std::cout<<"KLFitter::ReadConfigFile::ReadConfig(). Error: PathToInputFile==PathToOutputFile! Will not overwrite input file!!!"<<std::endl;return -1;}

}



bool KLFitter::ReadConfigFile::CheckIOPath()
{ 
  if(input_path==output_path)
    {return false;}
  else
    {return true;}
}

int KLFitter::ReadConfigFile::GetValue(double * ret, std::string line, size_t found)
{
  std::string::iterator k;
  bool is_whitespace=true;
  bool is_number=false;
  bool found_point=false;
  int a=0;
  int l=0;
  unsigned int i=0;
  *ret=0;
  if(line.size()<=found)
    {return -1;}
  else
    {

      for(k=line.begin();k<line.end();k++)
        {
          //go forward to the place of "="
          if(i<=found)
            {i++;continue;}
          //ignore any white space unless it is within a number
          if((*k==std::string(" ") || *k==std::string("\t")) && is_number==false)
            {is_whitespace=true;continue;}
          else
            {
              //Is this char a number (0 to 9) ?
              if(IsNumber(k,&a))
                {                                       
                  *ret+=a*pow(10,l);l--; // calculate the whole number. This result has to be corrected by a factor of some power of 10. But this correction is only known when the "." is found or the end of the number
                  is_number=true;
                                                
                }
              else
                {
                  if(*k==std::string(".") && found_point==false) //is this nonumber char the point?
                    {
                      *ret*=pow(10,fabs(l)-1); //if yes and there was no point before, correct the number and set the "l" this way, that it calculates the powers for the coming numbers in the right way
                      l=-1;     
                      found_point=true;                                                    
                    }
                  else
                    {
                      break;
                    }
                }
            }
        }
      if(found_point==false){*ret*=pow(10,fabs(l)-1);} //calculate the correction, if the number ended without a point
      if(is_number==true){return 1;}
      else{return -1;}
    }

}

int KLFitter::ReadConfigFile::GetPath(std::string * ret, std::string line, size_t found)
{
  std::string::iterator k;
  bool is_whitespace=true;
  bool path_started=false;
  int l=0;
  unsigned int i=0;
  *ret=std::string(""); //clear ret
  if(line.size()<=found)
    {return -1;}        
  else
    {

      for(k=line.begin();k<line.end();k++)
        {
          //go forward to the place of "="
          if(i<=found)
            {i++;continue;}
          //ignore any white space unless the path started
          if((*k==std::string(" ") || *k==std::string("\t")) && path_started==false)
            {is_whitespace=true;continue;}
          else
            {
              if(*k==std::string(" ") || *k==std::string("\t") || *k==std::string("#"))
                {break;}
              *ret+=*k;l++;
              path_started=true;
            }
        }
        
      if(l==0){*ret="";return -1;} //return -1 if path has length 0
      else{return 0;}           
    }

}

// --------------------------------------------------------- 
bool KLFitter::ReadConfigFile::IsNumber(std::string::iterator a, int * number)
{
        
  if(*a==std::string("0")){*number=0;return true;}
  else if (*a==std::string("1")){*number=1;return true;}
  else if (*a==std::string("2")){*number=2;return true;}
  else if (*a==std::string("3")){*number=3;return true;}
  else if (*a==std::string("4")){*number=4;return true;}
  else if (*a==std::string("5")){*number=5;return true;}
  else if (*a==std::string("6")){*number=6;return true;}
  else if (*a==std::string("7")){*number=7;return true;}
  else if (*a==std::string("8")){*number=8;return true;}
  else if (*a==std::string("9")){*number=9;return true;}
  else{return false;*number=-1;}
}


// --------------------------------------------------------- 
int KLFitter::ReadConfigFile::GetTrueOrFalse(std::string line, size_t found)
{
  if(line.size()<=found)
    {return -1;}
  else
    {
      //string to lower case

      for (unsigned int i = 0; i < line.size(); i++)
        {line[i] = tolower(line[i]);}   
      if(line.find("true", found)!=std::string::npos)
        {
          return 1;
        }
      else if(line.find("false", found)!=std::string::npos)
        {
          return 0;
        }
      else
        {return -1;}

    }
}

//<----------- End of reading out config.conf----------------->
