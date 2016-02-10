#include <sstream>

int drawx0(char* FILEN){

  TFile f( FILEN ) ;

  surfTuple->SetMarkerStyle( kOpenCircle ) ;
  surfTuple->SetMarkerColor( kBlue );
  matTuple->SetMarkerStyle( kPlus ) ;
  matTuple->SetMarkerColor( kRed );

  // get all theta and phi values
  std::vector<int> thetas, phis ;
  TH1F* ht=new TH1F("ht","ht",360,-0.5, 359.5 ) ;
  TH1F* hp=new TH1F("hp","hp",360,-0.5, 359.5 ) ;

  matTuple->Draw("theta>>ht") ;
  matTuple->Draw("phi>>hp") ;

  for(unsigned i=0,N=ht->GetNbinsX() ; i<N ; ++i ){
    if( ht->GetBinContent(i)  > 0.5 )  
      thetas.push_back(  int( ht->GetBinCenter(i) ) )  ; 
  }
  for(unsigned i=0,N=hp->GetNbinsX() ; i<N ; ++i ){
    if( hp->GetBinContent(i)  > 0.5 )  
      phis.push_back(  int( ht->GetBinCenter(i) ) )  ; 
  }

  std::string pdfFile( std::string( FILEN ) + std::string( ".pdf(" ) ) ;

  for(unsigned i=0,N=thetas.size();i<N;++i){
    for(unsigned j=0,M=phis.size();j<M;++j){

      if( i==N-1 && j==M-1 )
	pdfFile =  std::string( FILEN ) + std::string( ".pdf)" )  ;

      if( thetas[i] < 45 )
	draw_x0_z( thetas[i], phis[j] , pdfFile ) ;
      else
	draw_x0_r( thetas[i], phis[j] , pdfFile ) ;

    }
  }

}

int draw_x0_r(int theta, int phi , const std::string& pdfFile ){

  std::stringstream sCut ;
  sCut << "theta==" << theta << "&&phi==" << phi << "&&x0<.15" ;


  std::stringstream sTitle ;
  sTitle << "integrated x0 vs. r [theta=" << theta << ",phi=" << phi <<"]" ;
  
  c1 = new TCanvas("C1", sTitle.str().c_str()  ) ;
  
  c1->SetLogx(1) ;
  c1->SetLogy(0) ;
  
  // c1->Divide(1,2) ;
  // c1->cd(1) ;

  matTuple->Draw("x0:sqrt(epx*epx+epy*epy)", sCut.str().c_str() ) ;
  surfTuple->Draw("x0:sqrt(epx*epx+epy*epy)", sCut.str().c_str(),"same") ;


  //  std::string pdfFile( std::string( FILEN ) + std::string( ".pdf" ) ) ;

  std::cout << " writing to file " << pdfFile << " for theta : " << theta << " phi : "<< phi << std::endl ;
  c1->Print( pdfFile.c_str() ) ;
}


int draw_x0_z(int theta, int phi , const std::string& pdfFile ){

  std::stringstream sCut ;
  sCut << "theta==" << theta << "&&phi==" << phi << "&&x0<.15" ;


  std::stringstream sTitle ;
  sTitle << "integrated x0 vs. z [theta=" << theta << ",phi=" << phi <<"]" ;
  
  c1 = new TCanvas("C1", sTitle.str().c_str()  ) ;
  
  c1->SetLogx(1) ;
  c1->SetLogy(0) ;
  
  // c1->Divide(1,2) ;
  // c1->cd(1) ;

  matTuple->Draw( "x0:epz", sCut.str().c_str() ) ;
  surfTuple->Draw("x0:epz", sCut.str().c_str(),"same") ;


  //  std::string pdfFile( std::string( FILEN ) + std::string( ".pdf" ) ) ;

  std::cout << " writing to file " << pdfFile << " for theta : " << theta << " phi : "<< phi << std::endl ;
  c1->Print( pdfFile.c_str() ) ;
}
