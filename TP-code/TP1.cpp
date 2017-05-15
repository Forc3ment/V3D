#include <iostream>

#include <visp/vpDebug.h>
#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
#include <visp/vpImageSimulator.h>
#include <visp/vpDisplayX.h>


using namespace std ;



int main()
{
  vpImage<unsigned char> Ig(300,400,0);
  vpImage<unsigned char> Id(300,400,0);
  vpImage<vpRGBa> Iimage(876,1200);


  vpImageIo::read(Iimage,"../Data/big-sleep.jpg") ;


  double L = 0.600 ;
  double l = 0.438;
  // Initialise the 3D coordinates of the Iimage corners
  vpColVector X[4];
  for (int i = 0; i < 4; i++) X[i].resize(3);
  // Top left corner
  X[0][0] = -L;
  X[0][1] = -l;
  X[0][2] = 0;

  // Top right corner
  X[1][0] = L;
  X[1][1] = -l;
  X[1][2] = 0;

  // Bottom right corner
  X[2][0] = L;
  X[2][1] = l;
  X[2][2] = 0;

  //Bottom left corner
  X[3][0] = -L;
  X[3][1] = l;
  X[3][2] = 0;



  vpImageSimulator sim;
  sim.init(Iimage, X);
  vpCameraParameters cam(800.0, 800.0, 200, 150);




  cam.printParameters() ;


  // I1g
  vpHomogeneousMatrix  gMo(0,0,2,  vpMath::rad(0),vpMath::rad(0),0) ;
  sim.setCameraPosition(gMo);
  sim.getImage(Ig,cam);
  cout << "Image I1g " <<endl ;
  cout << "gMo " << endl ;
  cout << gMo << endl ;

  // I1d
  vpHomogeneousMatrix dMo(0.1,0,1.9,
			  vpMath::rad(5),vpMath::rad(5),vpMath::rad(5)) ;

  sim.setCameraPosition(dMo);
  sim.getImage(Id,cam);
  cout << "Image I1d " <<endl ;
  cout << "dMo " << endl ;
  cout << dMo << endl ;


  vpDisplayX dg(Ig,10,10,"Ig") ;
  vpDisplay::display(Ig) ;
  vpDisplay::flush(Ig) ;

  vpDisplayX dd(Id,10,10,"Id") ;
  vpDisplay::display(Id) ;
  vpDisplay::flush(Id) ;

  vpImageIo::write(Ig, "../result/Ig.pgm");
  vpImageIo::write(Id, "../result/Id.pgm");

  vpImagePoint pd, pg, p1, p2 ;

  double a, b, c, u1, u2, v1, v2;

  //Calcul de la matrice fondamentale

  //De droite vers gauche
  vpMatrix K_inverse = cam.get_K_inverse();
  vpMatrix K_inverse_transpose = K_inverse.transpose();

  vpHomogeneousMatrix gMd = gMo * dMo.inverse();
  vpTranslationVector gtd;
  gMd.extract(gtd);
  vpMatrix gTdx = gtd.skew();
  vpRotationMatrix gRd;
  gMd.extract(gRd);

  vpMatrix gFd = K_inverse_transpose * gTdx * gRd * K_inverse;

  //De gauche vers droite
  // vpMatrix K_inverse = cam.get_K_inverse();
  // vpMatrix K_inverse_transpose = K_inverse.transpose();
  //
  // vpHomogeneousMatrix dMg = dMo * gMo.inverse();
  // vpTranslationVector dtg;
  // dMg.extract(dtg);
  // vpMatrix dTgx = dtg.skew();
  // vpRotationMatrix dRg;
  // dMg.extract(dRg);
  //
  // vpMatrix dFg = K_inverse_transpose * dTgx * dRg * K_inverse;

  for (int i=0 ; i < 10 ; i++)
    {
      //De droite vers gauche
      cout << "Click point number " << i << endl ;
      vpDisplay::getClick(Id, pd) ;
      pd.set_u(50);
      pd.set_v(75);
      vpDisplay::displayCross(Id,pd,8,vpColor::green) ;


      // Calcul du lieu geometrique
      vpMatrix xd(3,1);
      xd[0][0] = pd.get_u();
      xd[1][0] = pd.get_v();
      xd[2][0] = 1;

      vpMatrix Deg = gFd * xd;

      a = Deg[0][0];
      b = Deg[1][0];
      c = Deg[2][0];

      u1 = 0;
      u2 = 399;
      std::cout << "V1 : " << v1 << " V2 : " << v2 << std::endl;
      std::cout << "a : " << a << " b : " << b << " c : " << c << std::endl;

      v1 = (-a*u1-c)/b;
      v2 = (-a*u2-c)/b;
      std::cout << "V1 : " << v1 << " V2 : " << v2 << std::endl;

      p1.set_u(u1);
      p1.set_v(v1);
      p2.set_u(u2);
      p2.set_v(v2);

      vpDisplay::displayLine(Ig, p1, p2, vpColor::green);

      vpDisplay::flush(Id) ;
      vpDisplay::flush(Ig) ;

      //De gauche vers droite
      // cout << "Click point number " << i << endl ;
      // vpDisplay::getClick(Ig, pg) ;
      // vpDisplay::displayCross(Ig,pg,8,vpColor::green) ;
      //
      //
      // // Calcul du lieu geometrique
      // vpMatrix xg(3,1);
      // xg[0][0] = pg.get_u();
      // xg[1][0] = pg.get_v();
      // xg[2][0] = 1;
      //
      // vpMatrix Deg = dFg * xg;
      //
      // a = Deg[0][0];
      // b = Deg[1][0];
      // c = Deg[2][0];
      //
      // u1 = 0;
      // u2 = 399;
      // std::cout << "V1 : " << v1 << " V2 : " << v2 << std::endl;
      // std::cout << "a : " << a << " b : " << b << " c : " << c << std::endl;
      //
      // v1 = (-a*u1-c)/b;
      // v2 = (-a*u2-c)/b;
      // std::cout << "V1 : " << v1 << " V2 : " << v2 << std::endl;
      //
      // p1.set_u(u1);
      // p1.set_v(v1);
      // p2.set_u(u2);
      // p2.set_v(v2);
      //
      // vpDisplay::displayLine(Id, p1, p2, vpColor::green);
      //
      // vpDisplay::flush(Id) ;
      // vpDisplay::flush(Ig) ;
    }

  // exemple de code pour sauvegarder une image avec les plan overlay
  vpImage<vpRGBa> Icol ;
  vpImage<vpRGBa> Ig_line ;
  vpDisplay::getImage(Id,Icol) ;
  vpDisplay::getImage(Ig,Ig_line) ;
  vpImageIo::write(Icol,"resultat.jpg") ;
  vpImageIo::write(Ig_line,"resultat_line.jpg") ;
  vpImageIo::write(Id,"I1g.jpg") ;




  vpDisplay::getClick(Id) ;
  cout << "OK " << endl ;

  vpDisplay::close(Id) ;
  vpDisplay::close(Ig) ;



  return 0;
}
