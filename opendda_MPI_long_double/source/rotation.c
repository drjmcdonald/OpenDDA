/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Functions to perform vector rotation

   Copyright (C) 2006 James Mc Donald
   Computational Astrophysics Laboratory
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "rotation.h"

void rotation_euler(long double phi,long double costheta,long double psi){
/* The negative ZYZ vector rotation matrix R(-psi)R(-theta)R(-phi)=
    _                    _
   | Rx_{x} Rx_{y} Rx_{z} |
   | Ry_{x} Ry_{y} Ry_{z} |
   |_Rz_{x} Rz_{y} Rz_{z}_|

   Rx_{x}=cos(phi)cos(theta)cos(psi)-sin(phi)sin(psi)
   Rx_{y}=sin(phi)cos(theta)cos(psi)+cos(phi)sin(psi)
   Rx_{z}=-sin(theta)cos(psi)
   Ry_{x}=-cos(phi)cos(theta)sin(psi)-sin(phi)cos(psi)
   Ry_{y}=-sin(phi)cos(theta)sin(psi)+cos(phi)cos(psi)
   Ry_{z}=sin(theta)sin(psi)
   Rz_{x}=cos(phi)sin(theta)
   Rz_{y}=sin(phi)sin(theta)
   Rz_{z}=cos(theta)

   The negative ZYZ vector rotation is given by:
          _                    _  _  _
         | Rx_{x} Rx_{y} Rx_{z} || xx |
         | Ry_{x} Ry_{y} Ry_{z} || xy |
         |_Rz_{x} Rz_{y} Rz_{z}_||_xz_| */

   long double cosphi,sinphi,theta,sintheta,cospsi,sinpsi,rotation[9];

   /* Euler theta is distributed in cos(theta) */
   theta=radians_to_degrees(acosl(costheta));
   cosphi=cosl(degrees_to_radians(phi)); /* Convert phi to radians and calculate the cosine */
   sinphi=sinl(degrees_to_radians(phi)); /* Convert phi to radians and calculate the sine */
   sintheta=sinl(degrees_to_radians(theta)); /* Convert theta to radians and calculate the sine */
   cospsi=cosl(degrees_to_radians(psi)); /* Convert psi to radians and calculate the cosine */
   sinpsi=sinl(degrees_to_radians(psi)); /* Convert psi to radians and calculate the sine */

   /* Perform a negative rotation */
   rotation[0]=cosphi*costheta*cospsi-sinphi*sinpsi; /* Rx_{x} */
   rotation[1]=sinphi*costheta*cospsi+cosphi*sinpsi; /* Rx_{y} */
   rotation[2]=-sintheta*cospsi; /* Rx_{z} */
   rotation[3]=-cosphi*costheta*sinpsi-sinphi*cospsi; /* Ry_{x} */
   rotation[4]=-sinphi*costheta*sinpsi+cosphi*cospsi; /* Ry_{y} */
   rotation[5]=sintheta*sinpsi; /* Ry_{z} */
   rotation[6]=cosphi*sintheta; /* Rz_{x} */
   rotation[7]=sinphi*sintheta; /* Rz_{y} */
   rotation[8]=costheta; /* Rz_{z} */

/* Rotate the incident vectors LF to TF */
   /* Rotate incident polarisation state LF to TF N.B. LF has no z component
   xx=incident_LF.polarisation[0]
   xy=incident_LF.polarisation[1]
   xz=incident_LF.polarisation[2] */
   incident_TF.polarisation[0]=ldcomplex_add(ldcomplex_scale(incident_LF.polarisation[0],rotation[0]),
                                 ldcomplex_scale(incident_LF.polarisation[1],rotation[1])); /* e0x in the target frame */
   incident_TF.polarisation[1]=ldcomplex_add(ldcomplex_scale(incident_LF.polarisation[0],rotation[3]),
                                 ldcomplex_scale(incident_LF.polarisation[1],rotation[4])); /* e0y in the target frame */
   incident_TF.polarisation[2]=ldcomplex_add(ldcomplex_scale(incident_LF.polarisation[0],rotation[6]),
                                 ldcomplex_scale(incident_LF.polarisation[1],rotation[7])); /* e0z in the target frame */
   /* Rotate orthonormal polarisation state LF to TF N.B. LF has no z component
   xx=incident_LF.orthonormal[0]
   xy=incident_LF.orthonormal[1]
   xz=incident_LF.orthonormal[2] */
   incident_TF.orthonormal[0]=ldcomplex_add(ldcomplex_scale(incident_LF.orthonormal[0],rotation[0]),
                              ldcomplex_scale(incident_LF.orthonormal[1],rotation[1])); /* e1x in the target frame */
   incident_TF.orthonormal[1]=ldcomplex_add(ldcomplex_scale(incident_LF.orthonormal[0],rotation[3]),
                              ldcomplex_scale(incident_LF.orthonormal[1],rotation[4])); /* e1y in the target frame */
   incident_TF.orthonormal[2]=ldcomplex_add(ldcomplex_scale(incident_LF.orthonormal[0],rotation[6]),
                              ldcomplex_scale(incident_LF.orthonormal[1],rotation[7])); /* e1z in the target frame */
   /* Rotate incident/propagation direction LF to TF N.B. incident_LF.n_inc=(0,0,1)
   xx=incident_LF.n_inc[0]
   xy=incident_LF.n_inc[1]
   xz=incident_LF.n_inc[2] */
   incident_TF.n_inc[0]=rotation[2]; /* n_{inc}x in the target frame */
   incident_TF.n_inc[1]=rotation[5]; /* n_{inc}y in the target frame */
   incident_TF.n_inc[2]=rotation[8]; /* n_{inc}z in the target frame */

/* Rotate the scattering vectors LF to TF */
   /* Rotate the scattering vector parallel to the scattering plane LF to TF
   xx=scattering_LF.theta_sca[0]
   xy=scattering_LF.theta_sca[1]
   xz=scattering_LF.theta_sca[2] */
   scattering_TF.theta_sca[0]=scattering_LF.theta_sca[0]*rotation[0]+scattering_LF.theta_sca[1]*rotation[1]+
                              scattering_LF.theta_sca[2]*rotation[2]; /* theta_{sca}x in the target frame */
   scattering_TF.theta_sca[1]=scattering_LF.theta_sca[0]*rotation[3]+scattering_LF.theta_sca[1]*rotation[4]+
                              scattering_LF.theta_sca[2]*rotation[5]; /* theta_{sca}y in the target frame */
   scattering_TF.theta_sca[2]=scattering_LF.theta_sca[0]*rotation[6]+scattering_LF.theta_sca[1]*rotation[7]+
                              scattering_LF.theta_sca[2]*rotation[8]; /* theta_{sca}z in the target frame */
   /* Rotate the scattering vector perpendicular to the scattering plane LF to TF
   xx=scattering_LF.phi_sca[0]
   xy=scattering_LF.phi_sca[1]
   xz=scattering_LF.phi_sca[2] */
   scattering_TF.phi_sca[0]=scattering_LF.phi_sca[0]*rotation[0]+scattering_LF.phi_sca[1]*rotation[1]+
                              scattering_LF.phi_sca[2]*rotation[2]; /* phi_{sca}x in the target frame */
   scattering_TF.phi_sca[1]=scattering_LF.phi_sca[0]*rotation[3]+scattering_LF.phi_sca[1]*rotation[4]+
                              scattering_LF.phi_sca[2]*rotation[5]; /* phi_{sca}y in the target frame */
   scattering_TF.phi_sca[2]=scattering_LF.phi_sca[0]*rotation[6]+scattering_LF.phi_sca[1]*rotation[7]+
                              scattering_LF.phi_sca[2]*rotation[8]; /* phi_{sca}z in the target frame */
   /* Rotate the scattering direction LF to TF
   xx=scattering_LF.n_sca[0]
   xy=scattering_LF.n_sca[1]
   xz=scattering_LF.n_sca[2] */
   scattering_TF.n_sca[0]=scattering_LF.n_sca[0]*rotation[0]+scattering_LF.n_sca[1]*rotation[1]+
                           scattering_LF.n_sca[2]*rotation[2]; /* n_{sca}x in the target frame */
   scattering_TF.n_sca[1]=scattering_LF.n_sca[0]*rotation[3]+scattering_LF.n_sca[1]*rotation[4]+
                           scattering_LF.n_sca[2]*rotation[5]; /* n_{sca}y in the target frame */
   scattering_TF.n_sca[2]=scattering_LF.n_sca[0]*rotation[6]+scattering_LF.n_sca[1]*rotation[7]+
                           scattering_LF.n_sca[2]*rotation[8]; /* n_{sca}z in the target frame */
}

void scattering_vectors(long double phi,long double theta){
/* Scattering directions using the spherical coodinates phi and theta
   (phi,theta)=minimum,increment,maximum

   The azimuthal angle phi is the clockwise angle [0,2pi], in the [x_{L}y_{L}]-plane,
   between the scattering plane and the +x_{L}-axis looking in the +z_{L}-axis direction
   i.e. the angle between the scattering plane and the xz-plane. The scattering plane is
   the plane containing the incident direction and the scattering direction [+z_{L}-axis
   & n_{sca}])
   phi=0 => scattering plane=xz-plane, phi=90 => scattering plane=yz-plane

   The polar/zenith angle theta is the angle [0,pi] in the scattering plane, between the
   incident direction +z_{L}-axis and the scattering direction n_{sca}
   theta=0 => forward scattering, scattering direction=+z_{L} axis
   theta=180 => backward scattering, scattering direction=-z_{L} axis

   Calculates the scattering vectors for the spherical coordinate angles phi and theta.
   i.e. Takes the lab frame and rotates it to coincide with the scattering direction

   x_{L} -> theta_{sca} i.e. vector parallel to the scattering plane in the lab frame
   y_{L} -> phi_{sca} i.e. vector perpendicular to the scattering plane in the lab frame
   z_{L} -> n_{sca} i.e. scattering direction in the lab frame

   The positive ZY axes rotation phi about z_{L} and theta about y^{phi} is given by
   matrix R(phi,theta) is R(theta)R(phi)
      _                    _
     | Rx_{x} Rx_{y} Rx_{z} |
     | Ry_{x} Ry_{y} Ry_{z} |
     |_Rz_{x} Rz_{y} Rz_{z}_|

     Rx_{x}=cos(phi)cos(theta)
     Rx_{y}=sin(phi)cos(theta)
     Rx_{z}=-sin(theta)
     Ry_{x}=sin(phi)
     Ry_{y}=cos(phi)
     Ry_{z}=0
     Rz_{x}=cos(phi)sin(theta)
     Rz_{y}=sin(phi)sin(theta)
     Rz_{z}=cos(theta)
 _                                    _                 _                    _
|theta_{sca}x theta_{sca}y theta_{sca}z|               | x_{L}x x_{L}y x_{L}z | 
| phi_{sca}x   phi_{sca}y   phi_{sca}z |=R(theta)R(phi)| y_{L}x y_{L}y y_{L}z |
|_ n_{sca}x     n_{sca}y     n_{sca}z _|               |_z_{L}x z_{L}y z_{L}z_|
                                                        _     _
                                                       | 1 0 0 |
                                        =R(theta)R(phi)| 0 1 0 |=R(theta)R(phi)I=R(theta)R(phi)
                                                       |_0 0 1_|
 _                                    _   _                                                 _
|theta_{sca}x theta_{sca}y theta_{sca}z| | cos(phi)cos(theta) sin(phi)cos(theta) -sin(theta) | 
| phi_{sca}x   phi_{sca}y   phi_{sca}z |=|    -sin(phi)           cos(phi)           0       |
|_ n_{sca}x     n_{sca}y     n_{sca}z _| |_cos(phi)sin(theta) sin(phi)sin(theta)  cos(theta)_| */

   long double cosphi,sinphi,costheta,sintheta;

   cosphi=cosl(degrees_to_radians(phi)); /* Convert phi to radians and calculate the cosine */
   sinphi=sinl(degrees_to_radians(phi)); /* Convert phi to radians and calculate the sine */
   costheta=cosl(degrees_to_radians(theta)); /* Convert theta to radians and calculate the cosine */
   sintheta=sinl(degrees_to_radians(theta)); /* Convert theta to radians and calculate the sine */

   scattering_LF.theta_sca[0]=cosphi*costheta; /* theta_{sca}x in the lab frame */
   scattering_LF.theta_sca[1]=sinphi*costheta; /* theta_{sca}y in the lab frame */
   scattering_LF.theta_sca[2]=-sintheta; /* theta_{sca}z in the lab frame */

   scattering_LF.phi_sca[0]=-sinphi; /* phi_{sca}x in the lab frame */
   scattering_LF.phi_sca[1]=cosphi; /* phi_{sca}y in the lab frame */
   scattering_LF.phi_sca[2]=0.0L; /* phi_{sca}z in the lab frame */

   scattering_LF.n_sca[0]=cosphi*sintheta; /* n_{sca}x in the lab frame */
   scattering_LF.n_sca[1]=sinphi*sintheta; /* n_{sca}y in the lab frame */
   scattering_LF.n_sca[2]=costheta; /* n_{sca}z in the lab frame */
}
