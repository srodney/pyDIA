#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))




void deconvolve3_columns(int width,int height,int rowstride,
                         float *data,float *buffer,float a,float b) {
    float *row;
    float q;
    int i, j;

    if (!height || !width)
        return;
	
    if (height == 1) {
        q = a + 2.0*b;
        for (j = 0; j < width; j++)
            data[j] /= q;
        return;
    }
    if (height == 2) {
        q = a*(a + 2.0*b);
        for (j = 0; j < width; j++) {
            buffer[0] = (a + b)/q*data[j] - b/q*data[rowstride + j];
            data[rowstride + j] = (a + b)/q*data[rowstride + j] - b/q*data[j];
            data[j] = buffer[0];
        }
        return;
    }
	
    /* Special-case first row */
    buffer[0] = a + b;
    /* Inner rows */
    for (i = 1; i < height-1; i++) {
        q = b/buffer[i-1];
        buffer[i] = a - q*b;
        row = data + (i - 1)*rowstride;
        for (j = 0; j < width; j++)
            row[rowstride + j] -= q*row[j];
    }
    /* Special-case last row */
    q = b/buffer[i-1];
    buffer[i] = a + b*(1.0 - q);
    row = data + (i - 1)*rowstride;
    for (j = 0; j < width; j++)
        row[rowstride + j] -= q*row[j];
    /* Go back */
    row += rowstride;
    for (j = 0; j < width; j++)
        row[j] /= buffer[i];
    do {
        i--;
        row = data + i*rowstride;
        for (j = 0; j < width; j++)
            row[j] = (row[j] - b*row[rowstride + j])/buffer[i];
    } while (i > 0);
}


void deconvolve3_rows(int width,int height,int rowstride,float *data,
                                 float *buffer,float a,float b) {
    float *row;
    float q;
    int i, j;

    if (!height || !width)
        return;
	
    if (width == 1) {
        q = a + 2.0*b;
        for (i = 0; i < height; i++)
            data[i*rowstride] /= q;
        return;
    }
    if (width == 2) {
        q = a*(a + 2.0*b);
        for (i = 0; i < height; i++) {
            row = data + i*rowstride;
            buffer[0] = (a + b)/q*row[0] - b/q*row[1];
            row[1] = (a + b)/q*row[1] - b/q*row[0];
            row[0] = buffer[0];
        }
        return;
    }
	
    /* Special-case first item */
    buffer[0] = a + b;
    /* Inner items */
    for (j = 1; j < width-1; j++) {
        q = b/buffer[j-1];
        buffer[j] = a - q*b;
        data[j] -= q*data[j-1];
    }
    /* Special-case last item */
    q = b/buffer[j-1];
    buffer[j] = a + b*(1.0 - q);
    data[j] -= q*data[j-1];
    /* Go back */
    data[j] /= buffer[j];
    do {
        j--;
        data[j] = (data[j] - b*data[j+1])/buffer[j];
    } while (j > 0);
	
    /* Remaining rows */
    for (i = 1; i < height; i++) {
        row = data + i*rowstride;
        /* Forward */
        for (j = 1; j < width-1; j++)
            row[j] -= b*row[j-1]/buffer[j-1];
        row[j] -= b*row[j-1]/buffer[j-1];
        /* Back */
        row[j] /= buffer[j];
        do {
            j--;
            row[j] = (row[j] - b*row[j+1])/buffer[j];
        } while (j > 0);
    }
}


void resolve_coeffs_2d(int width, int height, int rowstride, float *data) {
    float *buffer;
    int     max;
	
    max = width > height ? width : height;
    buffer = (float *)malloc(max*sizeof(float));
    deconvolve3_rows(width, height, rowstride, data, buffer, 13.0/21.0, 4.0/21.0);
    deconvolve3_columns(width, height, rowstride, data, buffer, 13.0/21.0, 4.0/21.0);
    free(buffer);
}


float interpolate_2d(float x,float y,int rowstride,float *coeff) {
    float wx[4], wy[4];
    int i, j;
    float v, vx;

    wx[0] = 4.0/21.0 + (-11.0/21.0 + (0.5 - x/6.0)*x)*x;
    wx[1] = 13.0/21.0 + (1.0/14.0 + (-1.0 + x/2.0)*x)*x;
    wx[2] = 4.0/21.0 + (3.0/7.0 + (0.5 - x/2.0)*x)*x;
    wx[3] = (1.0/42.0 + x*x/6.0)*x;
    wy[0] = 4.0/21.0 + (-11.0/21.0 + (0.5 - y/6.0)*y)*y;
    wy[1] = 13.0/21.0 + (1.0/14.0 + (-1.0 + y/2.0)*y)*y;
    wy[2] = 4.0/21.0 + (3.0/7.0 + (0.5 - y/2.0)*y)*y;
    wy[3] = (1.0/42.0 + y*y/6.0)*y;
	
    v = 0.0;
    for (i = 0; i < 4; i++) {
        vx = 0.0;
        for (j = 0; j < 4; j++)
            vx += coeff[i*rowstride + j]*wx[j];
        v += wy[i]*vx;
    }
	
    return v;
}


float integrated_profile(int profile_type, int idx, int idy, float xpos,
			 float ypos, float *psf_parameters, float *lut_0,
			 float *lut_xd, float *lut_yd) {

    int psf_size;
    float   psf_height, psf_sigma_x, psf_sigma_y, psf_xpos, psf_ypos;
    float   p0;
    int     ip, jp;
    float  pi=3.14159265,fwtosig=0.8493218;
    
    psf_size = (int) psf_parameters[0];
    psf_height = psf_parameters[1];
    psf_sigma_x = psf_parameters[2];
    psf_sigma_y = psf_parameters[3];
    psf_ypos = psf_parameters[4];
    psf_xpos = psf_parameters[5];

    if (profile_type == 0) {
    
       // gaussian

       // PSF at location (Idx,Idy). PSF is centred at (7.5,7.5)
       // Analytic part

       p0 = 0.5*psf_height*pi*fwtosig*fwtosig*
            (erf((idx-7.5+0.5)/(1.41421356*psf_sigma_x)) - 
             erf((idx-7.5-0.5)/(1.41421356*psf_sigma_x))) *
            (erf((idy-7.5+0.5)/(1.41421356*psf_sigma_y)) - 
             erf((idy-7.5-0.5)/(1.41421356*psf_sigma_y)));

       // Index into the lookup table

       ip = psf_size/2 + 2*idx - 15;
       jp = psf_size/2 + 2*idy - 15;
       if ((ip>=0) && (ip<=psf_size-1) && (jp>=0) && (jp<=psf_size-1)) {
          p0 += lut_0[ip+psf_size*jp] + lut_xd[ip+psf_size*jp]*(xpos-psf_xpos) +
                lut_yd[ip+psf_size*jp]*(ypos-psf_ypos);
       }

       return p0;

    } else if (profile_type == 1) {

       //  moffat25
       //  From iraf/noao/digiphot/daophot/daolib/profile.x

       float d[4][4] = {{ 0.0,         0.0,        0.0,        0.0},
                        {-0.28867513,  0.28867513, 0.0,        0.0},
                        {-0.38729833,  0.0,        0.38729833, 0.0},
                        {-0.43056816, -0.16999052, 0.16999052, 0.43056816}};
       float w[4][4] = {{1.0,         0.0,        0.0,        0.0},
                        {0.5,         0.5,        0.0,        0.0},
                        {0.27777778,  0.44444444, 0.27777778, 0.0},
                        {0.17392742,  0.32607258, 0.32607258, 0.17392742}};

       float alpha = 0.3195079;
       float  p1sq, p2sq, p1p2, dx, dy, xy, denom, func, x[4], xsq[4], p1xsq[4];
       float  y, ysq, p2ysq, wt, p4fod, wp4fod, wf;
       int    npt, ix, iy;

       p1sq = psf_parameters[2]*psf_parameters[2];
       p2sq = psf_parameters[3]*psf_parameters[3];
       p1p2 = psf_parameters[2]*psf_parameters[3];
       dx = idx-7.5+0.5;
       dy = idy-7.5+0.5;
       xy = dx * dy;
       
       denom = 1.0 + alpha * (dx*dx/p1sq + dy*dy/p2sq + xy*psf_parameters[4]);
       if (denom > 1.0e4) {
          return 0.0;
       }

       p0 = 0.0;
       func = 1.0 / (p1p2*pow(denom,2.5));
       if (func >= 0.046) {
          npt = 4;
       } else if (func >= 0.0022) {
          npt = 3;
       } else if (func >= 0.0001) {
          npt = 2;
       } else if (func >= 1.0e-10) {
          p0 = (2.5 - 1.0) * func;
       }

       if (func >= 0.0001) {
       
          for (ix=0; ix<npt; ix++) {
             x[ix] = dx + d[npt][ix];
             xsq[ix] = x[ix]*x[ix];
             p1xsq[ix] = xsq[ix]/p1sq;
          }

          for (iy=0; iy<npt; iy++) {
             y = dy + d[npt][iy];
             ysq = y*y;
             p2ysq = ysq/p2sq;
             for (ix=0; ix<npt; ix++) {
                wt = w[npt][iy] * w[npt][ix];
                xy = x[ix] * y;
                denom = 1.0 + alpha * (p1xsq[ix] + p2ysq + xy*psf_parameters[4]);
                func = (2.5 - 1.0) / (p1p2 * pow(denom,2.5) );
                p4fod = 2.5 * alpha * func / denom;
                wp4fod = wt * p4fod;
                wf = wt * func;
                p0 += wf;
             }
          }
          
       }

       p0 *= psf_parameters[1];

       // Index into the lookup table

       ip = psf_size/2 + 2*idx - 15;
       jp = psf_size/2 + 2*idy - 15;
       if ((ip>=0) && (ip<=psf_size-1) && (jp>=0) && (jp<=psf_size-1)) {
           p0 += lut_0[ip+psf_size*jp] + lut_xd[ip+psf_size*jp]*(xpos-psf_xpos) +
                 lut_yd[ip+psf_size*jp]*(ypos-psf_ypos);
       }
       
       return p0;

   } else {

      return 0.0;

   }
   
}




void cu_convolve_image_psf(int profile_type, int nx, int ny, int dx, int dy,
			   int dp, int ds, int n_coeff, int nkernel,
			   int kernel_radius,int *kxindex,
			   int *kyindex, int* ext_basis, float *psf_parameters,
			   float *psf_0, float *psf_xd, float *psf_yd, 
			   float *coeff,
			   float *cim1, float* cim2, float *tex0, float *tex1) {

   int     id, txa, tyb, txag, tybg,imax,jmax,idxmin,idxmax,idymin,idymax;
   int     np, ns, i, j, ii, ip, jp, ic, ki, a, b;
   int     d1, sidx, l, m, l1, m1, ig, jg;
   int     psf_size, ix, jx;
   float   x, y, p0, p1, p1g, cpsf_pixel;
   int     xpos, ypos;
   float   psf_height, psf_sigma_x, psf_sigma_y, psf_sigma_xy, psf_xpos, psf_ypos;
   float   gain,psf_rad,psf_rad2, px, py;
   float   sx2, sy2, sxy2, sx2msy2, sx2psy2; 
   float  psf_norm,dd;
   float  pi=3.14159265,fwtosig=0.8493218;

   float  psf_sum;
   float  cpsf[256];

   int     idx, idy, blockIdx, blockIdy, blockDimx=16, blockDimy=16, gridDimx, gridDimy;

   gridDimx = (nx-1)/dx+1;
   gridDimy = (ny-1)/dy+1;

   // number of polynomial coefficients per basis function
   np = (dp+1)*(dp+2)/2;
   ns = (ds+1)*(ds+2)/2;

   // PSF parameters
   psf_size = (int) psf_parameters[0];
   psf_height = psf_parameters[1];
   psf_sigma_x = psf_parameters[2];
   psf_sigma_y = psf_parameters[3];
   psf_ypos = psf_parameters[4];
   psf_xpos = psf_parameters[5];
   psf_rad = psf_parameters[6];
   gain = psf_parameters[7];
   if (psf_rad > 6.0) {
     psf_rad = 6.0;
   }
   psf_rad2 = psf_rad*psf_rad;


   // PSF integral
   psf_sum = 0.0;
   for (i=1; i<psf_size-1; i++) {
     for (j=1; j<psf_size-1; j++) {
       psf_sum += psf_0[i+j*psf_size];
     }
   }
 
   if (profile_type == 0) {
     // gaussian
     psf_norm = 0.25*psf_sum + psf_height*2*pi*fwtosig*fwtosig;
   } else if (profile_type == 1) {
     // moffat25
     psf_sigma_xy = psf_parameters[8];
     sx2 = psf_sigma_x*psf_sigma_x;
     sy2 = psf_sigma_y*psf_sigma_y;
     sxy2 = psf_sigma_xy*psf_sigma_xy;
     sx2msy2 = 1.0/sx2 - 1.0/sy2;
     sx2psy2 = 1.0/sx2 + 1.0/sy2;
     px = 1.0/sqrt( sx2psy2 + sqrt(sx2msy2*sx2msy2 + sxy2) );
     py = 1.0/sqrt( sx2psy2 - sqrt(sx2msy2*sx2msy2 + sxy2) );
     psf_norm = 0.25*psf_sum + psf_height*pi*(px*py)/(psf_sigma_x*psf_sigma_y);
   }
 


   for (blockIdx=0; blockIdx<gridDimx; blockIdx+=1) {
     for (blockIdy=0; blockIdy<gridDimy; blockIdy+=1) {

       // star position in normalised units
       xpos = blockIdx*dx + dx/2;
       ypos = blockIdy*dy + dy/2;
       x = (xpos - 0.5*(nx-1))/(nx-1);
       y = (ypos - 0.5*(ny-1))/(ny-1);


  
       // Construct the convolved PSF

       for (idx=0; idx<blockDimx; idx++) {
	 for (idy=0; idy<blockDimy; idy++) {
	   id = idx+idy*blockDimx;
	   cpsf[id] = 0.0;

	   // PSF at location (Idx,Idy). PSF is centred at (7.5,7.5)
	   // Analytic part

	   p0 = integrated_profile(profile_type, idx, idy, (float)xpos, (float)ypos,
				   psf_parameters, psf_0, psf_xd, psf_yd);


	   cpsf_pixel = 0.0;
    
	   // Iterate over coefficients

	   for (ic=0; ic<n_coeff; ic++) {
 
	     // basis function position
	     ki = ic < np ? 0 : (ic-np)/ns + 1;

	     if (ki<nkernel) {
       
	       a = kxindex[ki];
	       b = kyindex[ki];
 
	       // Set the polynomial degree for the subvector and the
	       // index within the subvector
	       if (ki == 0) {
		 d1 = dp;
		 sidx = ic;
	       } else {
		 d1 = ds;
		 sidx = ic - np - (ki-1)*ns;
	       }
       

	       // Compute the polynomial index (l,m) values corresponding
	       // to the index within the subvector
	       l1 = m1 = 0;
	       if (d1 > 0) {
		 i = 0;
		 for (l=0; l<=d1; l++) {
		   for (m=0; m<=d1-l; m++) {
		     if (i == sidx) {
		       l1 = l;
		       m1 = m;
		     }
		     i++;
		   }
		 }
	       }

	       // Indices into the PSF

	       if (ki > 0) {

		 txa = idx + a;
		 tyb = idy + b;

		 p1 = integrated_profile(profile_type, txa, tyb, (float)xpos, (float)ypos,
					 psf_parameters, psf_0, psf_xd, psf_yd);

		 // If we have an extended basis function, we need to
		 // average the PSF over a 3x3 grid
		 if (ext_basis[ki]) {

		   p1 = 0.0;
		   for (ig=-1; ig<2; ig++) {
		     for (jg=-1; jg<2; jg++) {
		       txag = txa + ig;
		       tybg = tyb + jg;
                               
		       p1g = integrated_profile(profile_type, txag, tybg, (float)xpos, 
						(float)ypos, psf_parameters, psf_0, 
						psf_xd, psf_yd);
		       p1 += p1g;
		     }
		   }
		   p1 /= 9.0;
           
		 }

		 cpsf_pixel += coeff[ic]*(p1-p0)*pow(x,l1)*pow(y,m1);

	       } else {

		 cpsf_pixel += coeff[ic]*p0*pow(x,l1)*pow(y,m1);

	       }
      
	     }
       
	   } //end ic loop


	   cpsf[id] = cpsf_pixel/psf_norm;

	 }
       }


       // Now convolve the image section with the convolved PSF

       imax = (xpos+dx/2) < nx ? xpos+dx/2 : nx;
       jmax = (ypos+dy/2) < ny ? ypos+dy/2 : ny;
       for (i=xpos-dx/2; i<imax; i++) {
	 for (j=ypos-dy/2; j<jmax; j++) {

	   cim1[i+j*nx] = 0.0;
	   cim2[i+j*nx] = 0.0;

	   idxmin = ((i-8) >= 0) ? 0 : 8-i;
	   idxmax = (blockDimx < (nx+8-i)) ? blockDimx : nx+8-i;
	   idymin = ((j-8) >= 0) ? 0 : 8-j;
	   idymax = (blockDimy < (ny+8-j)) ? blockDimx : ny+8-j;
	   for (idx=idxmin; idx<idxmax; idx++) {
	     for (idy=idymin; idy<idymax; idy++) {
	       id = idx+idy*blockDimx;
	       ix = i+idx-8;
	       jx = j+idy-8;
	       if ((ix < nx) && (jx < ny)) {
		 cim1[i+j*nx] += cpsf[id]*tex0[ix+nx*jx];
		 cim2[i+j*nx] += cpsf[id]*tex1[ix+nx*jx];
	       }

	     }
	   }

	 }
       }

     }
   }

   return;

}







void cu_photom(int profile_type,
	       int nx, int ny, int dp, int ds, int n_coeff, int nkernel,
	       int kernel_radius,int *kxindex,
	       int *kyindex, int* ext_basis, float *psf_parameters,
	       float *psf_0, float *psf_xd, float *psf_yd, float *posx,
	       float *posy, float *coeff, 
	       float *flux, float *dflux, long gridDimx, int blockDimx, 
	       int blockDimy,
	       float *tex0, float *tex1, int *group_boundaries, 
	       float *group_positions_x, float *group_positions_y, int ngroups) {

  int     id, txa, tyb, txag, tybg;
  int     np, ns, i, j, ip, jp, ic, ki, a, b;
  int     d1, sidx, l, m, l1, m1, ig, jg;
  int     psf_size, ix, jx;
  float   x, y, p0, p1, p1g, cpsf_pixel, xpos, ypos, dd;
  float   psf_height, psf_sigma_x, psf_sigma_y, psf_sigma_xy, psf_xpos, psf_ypos;
  float   psf_rad, psf_rad2, gain, fl, inv_var, px, py;
  float   sx2, sy2, sxy2, sx2msy2, sx2psy2; 
  float  subx, suby, psf_norm, bgnd;
  float  pi=3.14159265,fwtosig=0.8493218;

  float psf_sum;
  float cpsf[256],cpsf0[256],mpsf[256];
  float  fsum1, fsum2, fsum3;
  int idx, idy, i_group;
  long blockIdx, i_group_previous;

  printf("Doing photometry for %d groups\n",ngroups);

  // number of polynomial coefficients per basis function
  np = (dp+1)*(dp+2)/2;
  ns = (ds+1)*(ds+2)/2;

  // PSF parameters
  psf_size = (int) psf_parameters[0];
  psf_height = psf_parameters[1];
  psf_sigma_x = psf_parameters[2];
  psf_sigma_y = psf_parameters[3];
  psf_ypos = psf_parameters[4];
  psf_xpos = psf_parameters[5];
  psf_rad = psf_parameters[6];
  gain = psf_parameters[7];
  if (psf_rad > 6.0) {
    printf("Warning: resetting psf_rad to maximum value 6\n");
    psf_rad = 6.0;
  }
  psf_rad2 = psf_rad*psf_rad;

  // PSF integral  
  psf_sum = 0.0;
  for (i=1; i<psf_size-1; i++) {
    for (j=1; j<psf_size-1; j++) {
      psf_sum += psf_0[i+j*psf_size];
    }
  }

  if (profile_type == 0) {
    // gaussian
    psf_norm = 0.25*psf_sum + psf_height*2*pi*fwtosig*fwtosig;
  } else if (profile_type == 1) {
    // moffat25
    psf_sigma_xy = psf_parameters[8];
    sx2 = psf_sigma_x*psf_sigma_x;
    sy2 = psf_sigma_y*psf_sigma_y;
    sxy2 = psf_sigma_xy*psf_sigma_xy;
    sx2msy2 = 1.0/sx2 - 1.0/sy2;
    sx2psy2 = 1.0/sx2 + 1.0/sy2;
    px = 1.0/sqrt( sx2psy2 + sqrt(sx2msy2*sx2msy2 + sxy2) );
    py = 1.0/sqrt( sx2psy2 - sqrt(sx2msy2*sx2msy2 + sxy2) );
    psf_norm = 0.25*psf_sum + psf_height*pi*(px*py)/(psf_sigma_x*psf_sigma_y);
  }
   
  // Loop over star groups
  i_group_previous = 0;
  for (i_group=0; i_group<ngroups; i_group++) {

    // printf("processing star group %d\n",i_group);

    // Compute the PSF

    for (idx=0; idx<blockDimx; idx++) {
      for (idy=0; idy<blockDimy; idy++) {
	id = idx+idy*blockDimx;
	cpsf0[id] = 0.0;
      }
    }

    // star group position in normalised units
    xpos = group_positions_x[i_group];
    ypos = group_positions_y[i_group];
    x = (xpos - 0.5*(nx-1))/(nx-1);
    y = (ypos - 0.5*(ny-1))/(ny-1);


    // Construct the convolved PSF
    for (idx=0; idx<blockDimx; idx++) {
      for (idy=0; idy<blockDimy; idy++) {
	id = idx+idy*blockDimx;

   
	// PSF at location (Idx,Idy). PSF is centred at (7.5,7.5)
	// Analytic part

	p0 = integrated_profile(profile_type, idx, idy, xpos, ypos,
				psf_parameters, psf_0, psf_xd, psf_yd);

	cpsf_pixel = 0.0;
    

	// Convolve the PSF

	// Iterate over coefficients
	for (ic=0; ic<n_coeff; ic++) {
 
	  // basis function position
	  ki = ic < np ? 0 : (ic-np)/ns + 1;

	  if (ki<nkernel) {
       
	    a = kxindex[ki];
	    b = kyindex[ki];
 
	    // Set the polynomial degree for the subvector and the
	    // index within the subvector
	    if (ki == 0) {
	      d1 = dp;
	      sidx = ic;
	    } else {
	      d1 = ds;
	      sidx = ic - np - (ki-1)*ns;
	    }
       

	    // Compute the polynomial index (l,m) values corresponding
	    // to the index within the subvector
	    l1 = m1 = 0;
	    if (d1 > 0) {
	      i = 0;
	      for (l=0; l<=d1; l++) {
		for (m=0; m<=d1-l; m++) {
		  if (i == sidx) {
		    l1 = l;
		    m1 = m;
		  }
		  i++;
		}
	      }
	    }

	    // Indices into the PSF
	    if (ki > 0) {

	      txa = idx + a;
	      tyb = idy + b;

	      p1 = integrated_profile(profile_type, txa, tyb, xpos, ypos,
				      psf_parameters, psf_0, psf_xd, psf_yd);

	      // If we have an extended basis function, we need to
	      // average the PSF over a 3x3 grid
	      if (ext_basis[ki]) {

		p1 = 0.0;
		for (ig=-1; ig<2; ig++) {
		  for (jg=-1; jg<2; jg++) {
		    txag = txa + ig;
		    tybg = tyb + jg;
		    
		    p1g = integrated_profile(profile_type, txag, tybg, xpos, ypos,
					     psf_parameters, psf_0, psf_xd, 
					     psf_yd);
                        
		    p1 += p1g;
		  }
		}
		p1 /= 9.0;
		
	      }

	      cpsf_pixel += coeff[ic]*(p1-p0)*pow(x,l1)*pow(y,m1);

	    } else {

	      cpsf_pixel += coeff[ic]*p0*pow(x,l1)*pow(y,m1);

	    }
	    
	  } 
	  
	} //end ic loop

	cpsf0[id] = cpsf_pixel/psf_norm;

      }
    }

    printf("psf computed\n");

    // Loop over stars within the group
    for (blockIdx=i_group_previous; blockIdx<group_boundaries[i_group]; blockIdx++) { 

      // printf("processing star %d at (%8.2f,%8.2f)\n",blockIdx,posx[blockIdx],posy[blockIdx]);

      // Copy the PSF
      for (i=0; i<256; i++) cpsf[i] = cpsf0[i];

      // Convert PSF to cubic OMOMS representation
      resolve_coeffs_2d(16,16,16,cpsf);

      xpos = posx[blockIdx];
      ypos = posy[blockIdx];
      subx = ceil(xpos+0.5+0.0000000001) - (xpos+0.5);
      suby = ceil(ypos+0.5+0.0000000001) - (ypos+0.5);

      for (idx=0; idx<blockDimx; idx++) {
        for (idy=0; idy<blockDimy; idy++) {
          id = idx+idy*blockDimx;
          mpsf[id] = 0.0;
          if ((idx > 1) && (idx < 14) && (idy > 1) && (idy < 14)) {
            mpsf[id] = (float)interpolate_2d(subx,suby,16,&cpsf[idx-2+(idy-2)*blockDimx]);
          }
	}
      }

      fl = 0.0;
      for (j=0; j<3; j++) {
      
	fsum1 = fsum2 = fsum3 = 0.0;
      
	for (idx=0; idx<16; idx++) {
	  for (idy=0; idy<16; idy++) {
	    id = idx+idy*blockDimx;

	    if (pow(idx-7.5,2)+pow(idy-7.5,2) < psf_rad2) {


	      // Fit the mapped PSF to the difference image to compute an
	      // optimal flux estimate.
	      // Assume the difference image is in tex(:,:,0)
	      // and the inverse variance in tex(:,:,1).
	      // We have to iterate to get the correct variance.

	      ix = (int)floor(xpos+0.5)+idx-8.0;
	      jx = (int)floor(ypos+0.5)+idy-8.0;
	

	      inv_var = 1.0/(1.0/tex1[ix+nx*jx] + fl*mpsf[id]/gain);

	      fsum1 += mpsf[id]*tex0[ix+nx*jx]*inv_var;
	      fsum2 += mpsf[id]*mpsf[id]*inv_var;
	      fsum3 += mpsf[id];

	    }

      	  } 
	
	}
	fl = fsum1/fsum2;

      }
    
      // printf("flux for star %d: %f +/- %f\n",blockIdx,fl,sqrt(fsum3*fsum3/fsum2));
      flux[blockIdx] = fl;
      dflux[blockIdx] = sqrt(fsum3*fsum3/fsum2);

    }
    i_group_previous = group_boundaries[i_group];
  
  }

}


void cu_photom_converge(int profile_type,
			int nx, int ny, int dp, int ds, int n_coeff, int nkernel,
			int kernel_radius,int *kxindex,
			int *kyindex, int* ext_basis, float *psf_parameters,
			float *psf_0, float *psf_xd, float *psf_yd, float xpos,
			float ypos, float xpos0, 
			float ypos0, float *coeff, 
			float *flux, float *dflux,
			float *sjx1, float *sjx2, float *sjy1, float *sjy2, 
			int blockDimx, int blockDimy,
			float *tex0, float *tex1) {

  int     id, txa, tyb, txag, tybg;
  int     np, ns, i, j, ip, jp, ic, ki, a, b;
  int     d1, sidx, l, m, l1, m1, ig, jg;
  int     psf_size, ix, jx;
  float   x, y, p0, p1, p1g, cpsf_pixel, dd;
  float   psf_height, psf_sigma_x, psf_sigma_y, psf_sigma_xy, psf_xpos, psf_ypos;
  float  subx, suby, psf_norm, bgnd;

  float psf_sum;
  float cpsf[256],cpsf0[256];
  float  mpsf[256], psfxd[256], psfyd[256], psf_rad, psf_rad2;
  float  fsum1, fsum2, fsum3, fl, gain, inv_var;
  float   sx2, sy2, sxy2, sx2msy2, sx2psy2; 
  float  a1, sa1px, sa1py, spx, spy, spxy, px, py;
  int blockIdx, idx, idy, i_group, i_group_previous;
  float  pi=3.14159265,fwtosig=0.8493218;




  // number of polynomial coefficients per basis function
  np = (dp+1)*(dp+2)/2;
  ns = (ds+1)*(ds+2)/2;

  // PSF parameters
  psf_size = (int) psf_parameters[0];
  psf_height = psf_parameters[1];
  psf_sigma_x = psf_parameters[2];
  psf_sigma_y = psf_parameters[3];
  psf_ypos = psf_parameters[4];
  psf_xpos = psf_parameters[5];
  psf_rad = psf_parameters[6];
  gain = psf_parameters[7];
  if (psf_rad > 6.0) {
    printf("Warning: resetting psf_rad to maximum value 6\n");
    psf_rad = 6.0;
  }
  psf_rad2 = psf_rad*psf_rad;

  // PSF integral  
  psf_sum = 0.0;
  for (i=1; i<psf_size-1; i++) {
    for (j=1; j<psf_size-1; j++) {
      psf_sum += psf_0[i+j*psf_size];
    }
  }

  if (profile_type == 0) {
    // gaussian
    psf_norm = 0.25*psf_sum + psf_height*2*pi*fwtosig*fwtosig;
  } else if (profile_type == 1) {
    // moffat25
    psf_sigma_xy = psf_parameters[8];
    sx2 = psf_sigma_x*psf_sigma_x;
    sy2 = psf_sigma_y*psf_sigma_y;
    sxy2 = psf_sigma_xy*psf_sigma_xy;
    sx2msy2 = 1.0/sx2 - 1.0/sy2;
    sx2psy2 = 1.0/sx2 + 1.0/sy2;
    px = 1.0/sqrt( sx2psy2 + sqrt(sx2msy2*sx2msy2 + sxy2) );
    py = 1.0/sqrt( sx2psy2 - sqrt(sx2msy2*sx2msy2 + sxy2) );
    psf_norm = 0.25*psf_sum + psf_height*pi*(px*py)/(psf_sigma_x*psf_sigma_y);
  }


  // Compute the PSF

  for (idx=0; idx<blockDimx; idx++) {
    for (idy=0; idy<blockDimy; idy++) {
      id = idx+idy*blockDimx;
      cpsf0[id] = 0.0;
    }
  }

  // star position in normalised units
  x = (xpos0 - 0.5*(nx-1))/(nx-1);
  y = (ypos0 - 0.5*(ny-1))/(ny-1);

  // Construct the convolved PSF
  for (idx=0; idx<blockDimx; idx++) {
    for (idy=0; idy<blockDimy; idy++) {
      id = idx+idy*blockDimx;
      
   
      // PSF at location (Idx,Idy). PSF is centred at (7.5,7.5)
      // Analytic part

      p0 = integrated_profile(profile_type, idx, idy, xpos, ypos,
			      psf_parameters, psf_0, psf_xd, psf_yd);

      cpsf_pixel = 0.0;
    

      // Convolve the PSF

      // Iterate over coefficients
      for (ic=0; ic<n_coeff; ic++) {
 
	// basis function position
	ki = ic < np ? 0 : (ic-np)/ns + 1;
	
	if (ki<nkernel) {
       
	  a = kxindex[ki];
	  b = kyindex[ki];
 
	  // Set the polynomial degree for the subvector and the
	  // index within the subvector
	  if (ki == 0) {
	    d1 = dp;
	    sidx = ic;
	  } else {
	    d1 = ds;
	    sidx = ic - np - (ki-1)*ns;
	  }
       

	  // Compute the polynomial index (l,m) values corresponding
	  // to the index within the subvector
	  l1 = m1 = 0;
	  if (d1 > 0) {
	    i = 0;
	    for (l=0; l<=d1; l++) {
	      for (m=0; m<=d1-l; m++) {
		if (i == sidx) {
		  l1 = l;
		  m1 = m;
		}
		i++;
	      }
	    }
	  }

	  // Indices into the PSF
	  if (ki > 0) {
	    
	    txa = idx + a;
	    tyb = idy + b;

	    p1 = integrated_profile(profile_type, txa, tyb, xpos, ypos,
				    psf_parameters, psf_0, psf_xd, psf_yd);

	    // If we have an extended basis function, we need to
	    // average the PSF over a 3x3 grid
	    if (ext_basis[ki]) {

	      p1 = 0.0;
	      for (ig=-1; ig<2; ig++) {
		for (jg=-1; jg<2; jg++) {
		  txag = txa + ig;
		  tybg = tyb + jg;
		    
		  p1g = integrated_profile(profile_type, txag, tybg, xpos, ypos,
					   psf_parameters, psf_0, psf_xd, 
					   psf_yd);

		  p1 += p1g;
		}
	      }
	      p1 /= 9.0;
		
	    }

	    cpsf_pixel += coeff[ic]*(p1-p0)*pow(x,l1)*pow(y,m1);

	  } else {
	    
	    cpsf_pixel += coeff[ic]*p0*pow(x,l1)*pow(y,m1);

	  }
	    
	} 
	  
      } //end ic loop

      cpsf0[id] = cpsf_pixel/psf_norm;

    }
  }

    // Copy the PSF
  for (i=0; i<256; i++) cpsf[i] = cpsf0[i];


  // Convert PSF to cubic OMOMS representation
  resolve_coeffs_2d(16,16,16,cpsf);


  subx = ceil(xpos+0.5+0.0000000001) - (xpos+0.5);
  suby = ceil(ypos+0.5+0.0000000001) - (ypos+0.5);

  // Interpolate the PSF to the subpixel star coordinates
  for (i=0; i<256; i++) {
    mpsf[i] = 0.0;
    psfxd[i] = 0.0;
    psfyd[i] = 0.0;
  }

  psf_sum = -1.0e-6;
  for (idx=1; idx<15; idx++) {
    for (idy=1; idy<15; idy++) {
      id = idx+idy*blockDimx;
      mpsf[id] = (float)interpolate_2d(subx,suby,16,
				       &cpsf[(idx-2)+(idy-2)*blockDimx]);
      psfxd[id] = 0.0;
      psfyd[id] = 0.0;
      mpsf[id] = mpsf[id] > 0.0 ? mpsf[id] : 0.0;
      psf_sum += mpsf[id];
    }
  }

  // Compute PSF derivatives 
  for (idx=2; idx<14; idx++) {
    for (idy=2; idy<14; idy++) {
      id = idx+idy*blockDimx;
      psfxd[id] = 0.5*(mpsf[idx+1+idy*blockDimx]-mpsf[idx-1+idy*blockDimx]);
      psfyd[id] = 0.5*(mpsf[idx+(idy+1)*blockDimx]-mpsf[idx+(idy-1)*blockDimx]);
    }
  }


  if (psf_sum < 0.0) {
      *sjx1 = 0.0;
      *sjx2 = 1.0;
      *sjy1 = 0.0;
      *sjy2 = 1.0;
      *flux = 0.0;
      *dflux = 1.0e6;
  
  } else {

  // Compute optimal flux estimate
  fl = 0.0;
  for (j=0; j<3; j++) {
      
    fsum1 = fsum2 = fsum3 = 0.0;
    
    for (idx=0; idx<16; idx++) {
      for (idy=0; idy<16; idy++) {
	id = idx+idy*blockDimx;
	
	if (pow(idx-7.5,2)+pow(idy-7.5,2) < psf_rad2) {


	  // Fit the mapped PSF to the difference image to compute an
	  // optimal flux estimate.
	  // Assume the difference image is in tex(:,:,0)
	  // and the inverse variance in tex(:,:,1).
	  // We have to iterate to get the correct variance.

	  ix = (int)floor(xpos+0.5)+idx-8.0;
	  jx = (int)floor(ypos+0.5)+idy-8.0;

	  inv_var = 1.0/(1.0/tex1[ix+jx*32] + fl*mpsf[id]/gain);
	  if (isnan(inv_var)) {
	    printf("idx,idy,tex1,fl,mpsf: %d %d %f %f %f\n",idx,idy,tex1[ix+jx*32],fl,mpsf[id]);
	    for (idx=0; idx<16; idx++) {
	      for (idy=0; idy<16; idy++) {
		id = idx+idy*blockDimx;
		ix = (int)floor(xpos+0.5)+idx-8.0;
		jx = (int)floor(ypos+0.5)+idy-8.0;
		printf("idx,idy,tex1: %d %d %f\n",idx,idy,tex1[ix+jx*32]);
	      }
	    }
	    exit(1);
	  }

	  fsum1 += mpsf[id]*tex0[ix+jx*32]*inv_var;
	  fsum2 += mpsf[id]*mpsf[id]*inv_var;
	  fsum3 += mpsf[id];
	  
	}
	
      } 
      
    }
    fl = fsum1/fsum2;

    if (isnan(fl)) {
      printf("j: %d\n",j);
      for (idx=0; idx<16; idx++) {
	for (idy=0; idy<16; idy++) {
	  id = idx+idy*blockDimx;
	  if (pow(idx-7.5,2)+pow(idy-7.5,2) < psf_rad2) {
	    ix = (int)floor(xpos+0.5)+idx-8.0;
	    jx = (int)floor(ypos+0.5)+idy-8.0;
	    inv_var = 1.0/(1.0/tex1[ix+jx*32] + fl*mpsf[id]/gain);
	    printf("idx,idy,ix,jx,mpsf,inv_var,tex0: %d %d %d %d %f %f %f\n",idx,idy,ix,jx,mpsf[id],inv_var,tex0[ix+jx*32]);
	  }
	}
      }
      exit(1);
    }

  }
    
  *flux = fl;
  *dflux = sqrt(fsum3*fsum3/fsum2);

  // Compute contributions to coordinate correction 
  sa1px = sa1py = spx = spy = spxy = 0.0;
  for (idx=2; idx<14; idx++) {
    for (idy=2; idy<14; idy++) {
      id = idx+idy*blockDimx;
      ix = (int)floor(xpos+0.5)+idx-8.0;
      jx = (int)floor(ypos+0.5)+idy-8.0;
      inv_var = 1.0/(1.0/tex1[ix+jx*32] + fl*mpsf[id]/gain);
      a1 = tex0[ix+jx*32]-fl*mpsf[id];
      sa1px += a1*psfxd[id]*inv_var;
      sa1py += a1*psfyd[id]*inv_var;
      spx += psfxd[id]*psfxd[id]*inv_var;
      spy += psfyd[id]*psfyd[id]*inv_var;
      spxy += psfxd[id]*psfyd[id]*inv_var;
      if (isnan(a1)) {
	printf("tex0,fl,mpsf: %f %f %f\n",tex0[ix+jx*32],fl,mpsf[id]);
	printf("%d %d %f %f %f %f %f %f\n",idx,idy,a1,sa1px,sa1py,spx,spy,spxy);
       	exit(1);
      }
    }
  }
  *sjx1 = fl*(sa1px-sa1py*spxy/spy);
  *sjx2 = fl*fl*(spx-spxy*spxy/spy);
  *sjy1 = fl*(sa1py-sa1px*spxy/spx);
  *sjy2 = fl*fl*(spy-spxy*spxy/spx);
  //    printf("%f %f %f %f     %f %f\n",*sjx1,*sjx2,*sjy1,*sjy2,(*sjx1)/(*sjx2),(*sjy1)/(*sjy2));
  //  exit(1);       
  }
}




void cu_compute_model(int dp, int ds, int db, int *kxindex,
		      int *kyindex, int* ext_basis, int nkernel, float *coefficient,
		      float *M, int gridDimx, int gridDimy, 
		      float *tex0, float *tex1) {

  int  np, ns, nb, hs, idx, ki, a, b, d1, sidx, l, m, l1, m1, i;
  float x, y, Bi;
  int nx;

  int blockIdx, blockIdy;

  nx = gridDimx;

  // Calculate number of terms in subvectors
  np = (dp+1)*(dp+2)/2;
  ns = (ds+1)*(ds+2)/2;
  nb = (db+1)*(db+2)/2;
  hs = (nkernel-1)*ns+np+nb;

  for (blockIdx=0; blockIdx<gridDimx; blockIdx++) {
    for (blockIdy=0; blockIdy<gridDimy; blockIdy++) {

      x = (blockIdx - 0.5*(gridDimx-1))/(gridDimx-1);
      y = (blockIdy - 0.5*(gridDimy-1))/(gridDimy-1);
      
      for (idx = 0; idx < hs; idx++) {

	// This is the index of the subvector and its kernel offsets
	ki = idx < np ? 0 : (idx-np)/ns + 1;
	a = b = 0;
	if (ki<nkernel) {
	  a = kxindex[ki];
	  b = kyindex[ki];
	}

	if ((blockIdx+a>=0) && (blockIdx+a<gridDimx) && (blockIdy+b>=0) && (blockIdy+b<gridDimy)) {



	  // Set the polynomial degree for the subvector and the
	  // index within the subvector
	  if (ki == 0) {
	    d1 = dp;
	    sidx = idx;
	  } else if (ki < nkernel) {
	    d1 = ds;
	    sidx = idx - np - (ki-1)*ns;
	  } else {
	    d1 = db;
	    sidx = idx - np - (ki-1)*ns;
	  }

	  // Compute the (l,m) values corresponding to the index within
	  // the subvector
	  l1 = m1 = 0;
	  if (d1 > 0) {
	    i = 0;
	    for (l=0; l<=d1; l++) {
	      for (m=0; m<=d1-l; m++) {
		if (i == sidx) {
		  l1 = l;
		  m1 = m;
		}
		i++;
	      }
	    }
	  }

	  if (ki == 0) {
	    Bi = tex0[blockIdx+nx*blockIdy];
	  } else if (ki < nkernel) {
	    if (ext_basis[ki]) {
	      Bi = tex1[blockIdx+a+nx*(blockIdy+b)]-tex0[blockIdx+nx*blockIdy];
	    } else {
	      Bi = tex0[blockIdx+a+nx*(blockIdy+b)]-tex0[blockIdx+nx*blockIdy];
	    }
	  } else {
	    Bi = 1.0;
	  }
	
	  M[blockIdx+gridDimx*blockIdy] += coefficient[idx]*pow(x,l1)*pow(y,m1)*Bi;
	}
      }
    }
  }
}


void cu_compute_vector(int dp, int ds, int db, int nx,
		       int ny, int *kxindex, int *kyindex, int *ext_basis, int nkernel,
		       int kernelRadius,float *V, int BlockDimx, int gridDimx,
		       float *tex0, float *tex1, float *tex2, float *tex3, float *tex4) {

  int idx; 
  int np, ns, ki, a, b, d1, i, j;
  int l, m, l1, m1;
  float py, x, y, Bi;
  float temp;

  int blockIdx;


  // Calculate number of terms in subvectors
  np = (dp+1)*(dp+2)/2;
  ns = (ds+1)*(ds+2)/2;


  for (blockIdx=0; blockIdx<gridDimx; blockIdx++) {

    // This is the index of the subvector and its kernel offsets
    ki = blockIdx < np ? 0 : (blockIdx-np)/ns + 1;

    a = b = 0;
    if (ki<nkernel) {
      a = kxindex[ki];
      b = kyindex[ki];
    }
   
    // Set the polynomial degrees for the submatrix and the
    // indices within the submatrix
    if (ki == 0) {
      d1 = dp;
      idx = blockIdx;
    } else if (ki < nkernel) {
      d1 = ds;
      idx = blockIdx - np - (ki-1)*ns;
    } else {
      d1 = db;
      idx = blockIdx - np - (ki-1)*ns;
    }

    // Compute the (l,m) values corresponding to the index within
    // the subvector
    i = 0;
    for (l=0; l<=d1; l++) {
      for (m=0; m<=d1-l; m++) {
	if (i == idx) {
	  l1 = l;
	  m1 = m;
	}
	i++;
      }
    }   

    // Compute the contribution to V from each image location.
    // Use individual threads to sum over columns.
    // tex[:,:,0] is the reference image,
    // tex[:,:,1] is the blurred reference image,
    // tex[:,:,2] is the target image,
    // tex[:,:,3] is the inverse variance,
    // tex[:,:,4] is the mask.
    // Bi is the basis image value.

    Bi = 1.0;
    V[blockIdx] = 0.0;

    for (j=kernelRadius; j<ny-kernelRadius; j++) {
      y = (j - 0.5*(ny-1))/(ny-1);
      py = pow(y,m1);
      for (i=kernelRadius; i<nx-kernelRadius; i++) {
	x = (i - 0.5*(nx-1))/(nx-1);
	if (ki == 0) {
	  Bi = tex0[i+nx*j];
	} else if (ki < nkernel) {
	  if (ext_basis[ki]) {
	    Bi = tex1[i+a+nx*(j+b)]-tex0[i+nx*j];
	  } else {
	    Bi = tex0[i+a+nx*(j+b)]-tex0[i+nx*j];
	  }
	} else {
	  Bi = 1.0;
	}
	V[blockIdx] += pow(x,l1)*py*Bi*tex2[i+nx*j]*tex3[i+nx*j]*tex4[i+nx*j];
      }
    }

  }

}


void cu_compute_vector_stamps(int dp, int ds, int db, int nx, int ny, int nstamps, 
			      int stamp_half_width, float *stamp_xpos, float* stamp_ypos,
			      int *kxindex, int *kyindex, int *ext_basis, int nkernel,
			      int kernelRadius,float *V, int BlockDimx, int gridDimx,
			      float *tex0, float *tex1, float *tex2, float *tex3, float *tex4){

  int idx; 
  int np, ns, ki, a, b, d1, i, j, i1, i2, j1, j2;
  int l, m, l1, m1;
  float py, x, y, Bi;
  float temp;
   
  int blockIdx;

  // Calculate number of terms in subvectors
  np = (dp+1)*(dp+2)/2;
  ns = (ds+1)*(ds+2)/2;

  for (blockIdx=0; blockIdx<gridDimx; blockIdx++) {

    // This is the index of the subvector and its kernel offsets
    ki = blockIdx < np ? 0 : (blockIdx-np)/ns + 1;

    a = b = 0;
    if (ki<nkernel) {
      a = kxindex[ki];
      b = kyindex[ki];
    }
   
    // Set the polynomial degrees for the submatrix and the
    // indices within the submatrix
    if (ki == 0) {
      d1 = dp;
      idx = blockIdx;
    } else if (ki < nkernel) {
      d1 = ds;
      idx = blockIdx - np - (ki-1)*ns;
    } else {
      d1 = db;
      idx = blockIdx - np - (ki-1)*ns;
    }

    // Compute the (l,m) values corresponding to the index within
    // the subvector
    i = 0;
    for (l=0; l<=d1; l++) {
      for (m=0; m<=d1-l; m++) {
	if (i == idx) {
	  l1 = l;
	  m1 = m;
	}
	i++;
      }
    }   

    // Compute the contribution to V from each image location.
    // Use individual threads to sum over columns.
    // tex[:,:,0] is the reference image,
    // tex[:,:,1] is the blurred reference image,
    // tex[:,:,2] is the target image,
    // tex[:,:,3] is the inverse variance,
    // tex[:,:,4] is the mask.
    // Bi is the basis image value.

    Bi = 1.0;
    V[blockIdx] = 0.0;

    for (idx = 0; idx<nstamps; idx++) {
      j1 = max(0,(int)stamp_ypos[idx]-stamp_half_width);
      j2 = min(ny,(int)stamp_ypos[idx]+stamp_half_width);
      for (j=j1; j<j2; j++) {
	y = (j - 0.5*(ny-1))/(ny-1);
	py = pow(y,m1);
	i1 = max(0,(int)stamp_xpos[idx]-stamp_half_width);
	i2 = min(nx,(int)stamp_xpos[idx]+stamp_half_width);
	for (i=i1; i<i2; i++) {
	  x = (i - 0.5*(nx-1))/(nx-1);
	  if (ki == 0) {
	    Bi = tex0[i+nx*j];
	  } else if (ki < nkernel) {
	    if (ext_basis[ki]) {
	      Bi = tex1[i+a+nx*(j+b)]-tex0[i+nx*j];
	    } else {
	      Bi = tex0[i+a+nx*(j+b)]-tex0[i+nx*j];
	    }
	  } else {
	    Bi = 1.0;
	  }
	  V[blockIdx] += pow(x,l1)*py*Bi*tex2[i+nx*j]*tex3[i+nx*j]*tex4[i+nx*j];
	}
      }
    }

  }
}


void cu_compute_matrix(int dp, int ds, int db, int nx, int ny, int *kxindex, 
		       int *kyindex, int *ext_basis, int nkernel, int kernelRadius,
		       float *H, int BlockDimx, int gridDimx, int gridDimy,
		       float *tex0, float *tex1, float *tex3, float *tex4){

  int idx, idy, idx0, idy0, idx1, idy1; 
  int np, ns, ki, kj, a, b, c, d, d1, d2, i, j;
  int l, m, l1, m1, l2, m2;
  float *px, *py, *x, *y, Bi, Bj,temp, pyj;
  float *ppx;
  float  *pBi1, *pBi2, *pBj1, *pBj2, *pt3, *pt4;


  int blockIdx, blockIdy;


  x = (float *) malloc(nx * sizeof(float));
  y = (float *) malloc(ny * sizeof(float));
  px = (float *) malloc(nx * sizeof(float));
  py = (float *) malloc(ny * sizeof(float));

  for (j=0; j<nx; j++) {
    x[j] = (j - 0.5*(nx-1))/(nx-1);
  }
  for (j=0; j<ny; j++) {
    y[j] = (j - 0.5*(ny-1))/(ny-1);
  }

  for (blockIdy=0; blockIdy<gridDimy; blockIdy++) {
    for (blockIdx=0; blockIdx<gridDimx; blockIdx++) {
      H[blockIdx+gridDimx*blockIdy] = 0.0;
    }
  }

  // Calculate number of terms in submatrices
  np = (dp+1)*(dp+2)/2;
  ns = (ds+1)*(ds+2)/2;


  for (blockIdy=0; blockIdy<gridDimy; blockIdy++) {
    for (blockIdx=0; blockIdx<=blockIdy; blockIdx++) {


      // These are indices of the submatrix and their kernel offsets
      ki = blockIdx < np ? 0 : (blockIdx-np)/ns + 1;
      kj = blockIdy < np ? 0 : (blockIdy-np)/ns + 1;

      a = b = 0;
      if (ki<nkernel) {
	a = kxindex[ki];
	b = kyindex[ki];
      }
      if (kj<nkernel) {
	c = kxindex[kj];
	d = kyindex[kj];
      }
   
      // Set the polynomial degrees for the submatrix and the
      // indices within the submatrix
      if (ki == 0) {
	d1 = dp;
	idx = blockIdx;
      } else if (ki < nkernel) {
	d1 = ds;
	idx = blockIdx - np - (ki-1)*ns;
      } else {
	d1 = db;
	idx = blockIdx - np - (ki-1)*ns;
      }
      if (kj == 0) {
	d2 = dp;
	idy = blockIdy;
      } else if (kj < nkernel) {
	d2 = ds;
	idy = blockIdy - np - (kj-1)*ns;
      } else {
	d2 = db;
	idy = blockIdy - np - (kj-1)*ns;
      }

      if ((ki>0) && (ki<nkernel) && (kj>0) && (kj<nkernel) && (idx > idy)) {
	continue;
      }

      idx0 = idx;
      idy0 = idy;
      
      // Compute the (l,m) values corresponding to the indices within
      // the submatrix
      i = 0;
      for (l=0; l<=d1; l++) {
	for (m=0; m<=d1-l; m++) {
	  if (i == idx) {
	    l1 = l;
	    m1 = m;
	  }
	  i++;
	}
      }   
      i = 0;
      for (l=0; l<=d2; l++) {
	for (m=0; m<=d2-l; m++) {
	  if (i == idy) {
	    l2 = l;
	    m2 = m;
	  }
	  i++;
	}
      }

      for (j=0; j<nx; j++) {
	px[j] = pow(x[j],l1+l2);
      }
      for (j=0; j<ny; j++) {
	py[j] = pow(y[j],m1+m2);
      }

      // Compute the contribution to H from each image location.
      // Use individual threads to sum over columns.
      // tex[:,:,0] is the reference image,
      // tex[:,:,1] is the blurred reference image,
      // tex[:,:,2] is the target image,
      // tex[:,:,3] is the inverse variance,
      // tex[:,:,4] is the mask.
      // Bi and Bj are the basis image values.

      Bi = Bj = 1.0;
      temp =0.0;
      

      // # The following lines are an unravelling of this code with pointers
      //
      //      for (j=kernelRadius; j<ny-kernelRadius; j++) {
      //	for (i=kernelRadius; i<nx-kernelRadius; i++) {
      //	  if (ki == 0) {
      //	    Bi = tex0[i+nx*j];
      //	  } else if (ki < nkernel) {
      //	    if (ext_basis[ki]) {
      //	      Bi = tex1[i+a+nx*(j+b)]-tex0[i+nx*j];
      //	    } else {
      //	      Bi = tex0[i+a+nx*(j+b)]-tex0[i+nx*j];
      //	    }
      //	  } else {
      //	    Bi = 1.0;
      //	  }
      //	  if (kj == 0) {
      //	    Bj = tex0[i+nx*j];
      //	  } else if (kj < nkernel) {
      //	    if (ext_basis[kj]) {
      //	      Bj = tex1[i+c+nx*(j+d)]-tex0[i+nx*j];
      //	    } else {
      //	      Bj = tex0[i+c+nx*(j+d)]-tex0[i+nx*j];
      //	    }
      //	  } else {
      //	    Bj = 1.0;
      //	  }
      //	  temp += px[i]*py[j]*Bi*Bj*tex3[i+nx*j]*tex4[i+nx*j];
      //	}
      //      }



      if (ki == 0) {

	if (kj == 0) {

	  for (j=kernelRadius; j<ny-kernelRadius; j++) {
	    pBi1 = &tex0[kernelRadius+nx*j];
	    pBj1 = &tex0[kernelRadius+nx*j];
	    pt3 = &tex3[kernelRadius+nx*j];
	    pt4 = &tex4[kernelRadius+nx*j];
	    ppx = &px[kernelRadius];
	    pyj = py[j];
	    for (i=kernelRadius; i<nx-kernelRadius; i++) {
	      temp += (*ppx++)*pyj*(*pBi1++)*(*pBj1++)*(*pt3++)*(*pt4++);
	    }
	  }

	} else if (kj < nkernel) {

	  for (j=kernelRadius; j<ny-kernelRadius; j++) {
	    pBi1 = &tex0[kernelRadius+nx*j];
	    pBj1 = &tex0[kernelRadius+nx*j];
	    if (ext_basis[kj]) {
	      pBj2 = &tex1[kernelRadius+c+nx*(j+d)];
	    } else {
	      pBj2 = &tex0[kernelRadius+c+nx*(j+d)];
	    }
	    pt3 = &tex3[kernelRadius+nx*j];
	    pt4 = &tex4[kernelRadius+nx*j];
	    ppx = &px[kernelRadius];
	    pyj = py[j];
	    for (i=kernelRadius; i<nx-kernelRadius; i++) {
	      temp += (*ppx++)*pyj*(*pBi1++)*((*pBj2++)-(*pBj1++))*(*pt3++)*(*pt4++);
	    }
	  }

	} else {

	  for (j=kernelRadius; j<ny-kernelRadius; j++) {
	    pBi1 = &tex0[kernelRadius+nx*j];
	    pt3 = &tex3[kernelRadius+nx*j];
	    pt4 = &tex4[kernelRadius+nx*j];
	    ppx = &px[kernelRadius];
	    pyj = py[j];
	    for (i=kernelRadius; i<nx-kernelRadius; i++) {
	      temp += (*ppx++)*pyj*(*pBi1++)*(*pt3++)*(*pt4++);
	    }
	  }

	} 


      } else if (ki < nkernel) {

	if (kj == 0) {

	  for (j=kernelRadius; j<ny-kernelRadius; j++) {
	    pBi1 = &tex0[kernelRadius+nx*j];
	    if (ext_basis[ki]) {
	      pBi2 = &tex1[kernelRadius+a+nx*(j+b)];
	    } else {
	      pBi2 = &tex0[kernelRadius+a+nx*(j+b)];
	    }
	    pBj1 = &tex0[kernelRadius+nx*j];
	    pt3 = &tex3[kernelRadius+nx*j];
	    pt4 = &tex4[kernelRadius+nx*j];
	    ppx = &px[kernelRadius];
	    pyj = py[j];
	    for (i=kernelRadius; i<nx-kernelRadius; i++) {
	      temp += (*ppx++)*pyj*((*pBi2++)-(*pBi1++))*(*pBj1++)*(*pt3++)*(*pt4++);
	    }
	  }

	} else if (kj < nkernel) {

	  for (j=kernelRadius; j<ny-kernelRadius; j++) {
	    pBi1 = &tex0[kernelRadius+nx*j];
	    pBj1 = &tex0[kernelRadius+nx*j];
	    if (ext_basis[ki]) {
	      pBi2 = &tex1[kernelRadius+a+nx*(j+b)];
	    } else {
	      pBi2 = &tex0[kernelRadius+a+nx*(j+b)];
	    }
	    if (ext_basis[kj]) {
	      pBj2 = &tex1[kernelRadius+c+nx*(j+d)];
	    } else {
	      pBj2 = &tex0[kernelRadius+c+nx*(j+d)];
	    }
	    pt3 = &tex3[kernelRadius+nx*j];
	    pt4 = &tex4[kernelRadius+nx*j];
	    ppx = &px[kernelRadius];
	    pyj = py[j];
	    for (i=kernelRadius; i<nx-kernelRadius; i++) {
	      temp += (*ppx++)*pyj*((*pBi2++)-(*pBi1++))*((*pBj2++)-(*pBj1++))*(*pt3++)*(*pt4++);
	    }
	  }

	} else {

	  for (j=kernelRadius; j<ny-kernelRadius; j++) {
	    pBi1 = &tex0[kernelRadius+nx*j];
	    if (ext_basis[ki]) {
	      pBi2 = &tex1[kernelRadius+a+nx*(j+b)];
	    } else {
	      pBi2 = &tex0[kernelRadius+a+nx*(j+b)];
	    }
	    pt3 = &tex3[kernelRadius+nx*j];
	    pt4 = &tex4[kernelRadius+nx*j];
	    ppx = &px[kernelRadius];
	    pyj = py[j];
	    for (i=kernelRadius; i<nx-kernelRadius; i++) {
	      temp += (*ppx++)*pyj*((*pBi2++)-(*pBi1++))*(*pt3++)*(*pt4++);
	    }
	  }

	}


      } else {

	if (kj == 0) {

	  for (j=kernelRadius; j<ny-kernelRadius; j++) {
	    pBj1 = &tex0[kernelRadius+nx*j];
	    pt3 = &tex3[kernelRadius+nx*j];
	    pt4 = &tex4[kernelRadius+nx*j];
	    ppx = &px[kernelRadius];
	    pyj = py[j];
	    for (i=kernelRadius; i<nx-kernelRadius; i++) {
	      temp += (*ppx++)*pyj*(*pBj1++)*(*pt3++)*(*pt4++);
	    }
	  }

	} else if (kj < nkernel) {

	  for (j=kernelRadius; j<ny-kernelRadius; j++) {
	    pBj1 = &tex0[kernelRadius+nx*j];
	    if (ext_basis[kj]) {
	      pBj2 = &tex1[kernelRadius+c+nx*(j+d)];
	    } else {
	      pBj2 = &tex0[kernelRadius+c+nx*(j+d)];
	    }
	    pt3 = &tex3[kernelRadius+nx*j];
	    pt4 = &tex4[kernelRadius+nx*j];
	    ppx = &px[kernelRadius];
	    pyj = py[j];
	    for (i=kernelRadius; i<nx-kernelRadius; i++) {
	      temp += (*ppx++)*pyj*((*pBj2++)-(*pBj1++))*(*pt3++)*(*pt4++);
	    }
	  }

	} else {

	  for (j=kernelRadius; j<ny-kernelRadius; j++) {
	    pt3 = &tex3[kernelRadius+nx*j];
	    pt4 = &tex4[kernelRadius+nx*j];
	    ppx = &px[kernelRadius];
	    pyj = py[j];
	    for (i=kernelRadius; i<nx-kernelRadius; i++) {
	      temp += (*ppx++)*pyj*(*pt3++)*(*pt4++);
	    }
	  }

	}

      }

      H[blockIdx+gridDimx*blockIdy] = temp;
      H[blockIdy+gridDimx*blockIdx] = temp;
      if ((ki>0) && (ki<nkernel) && (kj>0) && (kj<nkernel)) {
	idx1 = np + (ki-1)*ns;
	idy1 = np + (kj-1)*ns;
	H[(idx1+idy0)+gridDimx*(idy1+idx0)] = temp;
	H[(idy1+idx0)+gridDimx*(idx1+idy0)] = temp;
      }

    }
  }

}



void cu_compute_matrix_stamps(int dp, int ds, int db, int nx, int ny, int nstamps, 
			      int stamp_half_width, float *stamp_xpos, float* stamp_ypos,
			      int *kxindex, int *kyindex, int *ext_basis, int nkernel,
			      int kernelRadius,float *H, int BlockDimx, int gridDimx, 
			      int gridDimy,
			      float *tex0, float *tex1, float *tex3, float *tex4) {

  int idx, idy, idx0, idy0, idx1, idy1; 
  int np, ns, ki, kj, a, b, c, d, d1, d2, i, j, i1, i2, j1, j2;
  int l, m, l1, m1, l2, m2;
  float *px, *py, *x, *y, Bi, Bj,temp, pyj;
  float *ppx;
  float  *pBi1, *pBi2, *pBj1, *pBj2, *pt3, *pt4;

  int blockIdx, blockIdy;

  x = (float *) malloc(nx * sizeof(float));
  y = (float *) malloc(ny * sizeof(float));
  px = (float *) malloc(nx * sizeof(float));
  py = (float *) malloc(ny * sizeof(float));

  for (j=0; j<nx; j++) {
    x[j] = (j - 0.5*(nx-1))/(nx-1);
  }
  for (j=0; j<ny; j++) {
    y[j] = (j - 0.5*(ny-1))/(ny-1);
  }

  for (blockIdy=0; blockIdy<gridDimy; blockIdy++) {
    for (blockIdx=0; blockIdx<gridDimx; blockIdx++) {
      H[blockIdx+gridDimx*blockIdy] = 0.0;
    }
  }


  // Calculate number of terms in submatrices
  np = (dp+1)*(dp+2)/2;
  ns = (ds+1)*(ds+2)/2;


  for (blockIdy=0; blockIdy<gridDimy; blockIdy++) {
    for (blockIdx=0; blockIdx<=blockIdy; blockIdx++) {

      // These are indices of the submatrix and their kernel offsets
      ki = blockIdx < np ? 0 : (blockIdx-np)/ns + 1;
      kj = blockIdy < np ? 0 : (blockIdy-np)/ns + 1;

      a = b = 0;
      if (ki<nkernel) {
	a = kxindex[ki];
	b = kyindex[ki];
      }
      if (kj<nkernel) {
	c = kxindex[kj];
	d = kyindex[kj];
      }
   
      // Set the polynomial degrees for the submatrix and the
      // indices within the submatrix
      if (ki == 0) {
	d1 = dp;
	idx = blockIdx;
      } else if (ki < nkernel) {
	d1 = ds;
	idx = blockIdx - np - (ki-1)*ns;
      } else {
	d1 = db;
	idx = blockIdx - np - (ki-1)*ns;
      }
      if (kj == 0) {
	d2 = dp;
	idy = blockIdy;
      } else if (kj < nkernel) {
	d2 = ds;
	idy = blockIdy - np - (kj-1)*ns;
      } else {
	d2 = db;
	idy = blockIdy - np - (kj-1)*ns;
      }

      if ((ki>0) && (ki<nkernel) && (kj>0) && (kj<nkernel) && (idx > idy)) {
	continue;
      }

      idx0 = idx;
      idy0 = idy;

      // Compute the (l,m) values corresponding to the indices within
      // the submatrix
      i = 0;
      for (l=0; l<=d1; l++) {
	for (m=0; m<=d1-l; m++) {
	  if (i == idx) {
	    l1 = l;
	    m1 = m;
	  }
	  i++;
	}
      }   
      i = 0;
      for (l=0; l<=d2; l++) {
	for (m=0; m<=d2-l; m++) {
	  if (i == idy) {
	    l2 = l;
	    m2 = m;
	  }
	  i++;
	}
      }

      for (j=0; j<nx; j++) {
	px[j] = pow(x[j],l1+l2);
      }
      for (j=0; j<ny; j++) {
	py[j] = pow(y[j],m1+m2);
      }

      // Compute the contribution to H from each image location.
      // Use individual threads to sum over stamps.
      // tex[:,:,0] is the reference image,
      // tex[:,:,1] is the blurred reference image,
      // tex[:,:,2] is the target image,
      // tex[:,:,3] is the inverse variance,
      // tex[:,:,4] is the mask.
      // Bi and Bj are the basis image values.

      Bi = Bj = 1.0;

      temp = 0.0;

      // # The following lines are an unravelling of this code with pointers

      //for (idx = 0; idx<nstamps; idx++) {
      //i1 = max(0,(int)stamp_xpos[idx]-stamp_half_width);
      //i2 = min(nx,(int)stamp_xpos[idx]+stamp_half_width);
      //for (i=i1; i<i2; i++) {
      //  x = (i - 0.5*(nx-1))/(nx-1);
      //  px = pow(x,l1+l2);
      //  j1 = max(0,(int)stamp_ypos[idx]-stamp_half_width);
      //  j2 = min(ny,(int)stamp_ypos[idx]+stamp_half_width);
      //  for (j=j1; j<j2; j++) {
      //    y = (j - 0.5*(ny-1))/(ny-1);
      //    py = pow(y,m1+m2);
      //    if (ki == 0) {
      //      Bi = tex0[i+nx*j];
      //    } else if (ki < nkernel) {
      //      if (ext_basis[ki]) {
      //	Bi = tex1[i+a+nx*(j+b)]-tex0[i+nx*j];
      //      } else {
      //	Bi = tex0[i+a+nx*(j+b)]-tex0[i+nx*j];
      //      }
      //    } else {
      //      Bi = 1.0;
      //    }
      //    if (kj == 0) {
      //      Bj = tex0[i+nx*j];
      //    } else if (kj < nkernel) {
      //      if (ext_basis[kj]) {
      //	Bj = tex1[i+c+nx*(j+d)]-tex0[i+nx*j];
      //      } else {
      //	Bj = tex0[i+c+nx*(j+d)]-tex0[i+nx*j];
      //      }
      //    } else {
      //      Bj = 1.0;
      //    }
      //    temp += px*py*Bi*Bj*tex3[i+nx*j]*tex4[i+nx*j];
      //  }
      //}
      //}


      if (ki == 0) {

	if (kj == 0) {

	  for (idx = 0; idx<nstamps; idx++) {

	    j1 = max(0,(int)stamp_ypos[idx]-stamp_half_width);
	    j2 = min(ny,(int)stamp_ypos[idx]+stamp_half_width);
	    i1 = max(0,(int)stamp_xpos[idx]-stamp_half_width);
	    i2 = min(nx,(int)stamp_xpos[idx]+stamp_half_width);

	    for (j=j1; j<j2; j++) {
	      pBi1 = &tex0[i1+nx*j];
	      pBj1 = &tex0[i1+nx*j];
	      pt3 = &tex3[i1+nx*j];
	      pt4 = &tex4[i1+nx*j];
	      ppx = &px[i1];
	      pyj = py[j];
	      for (i=i1; i<i2; i++) {
		temp += (*ppx++)*pyj*(*pBi1++)*(*pBj1++)*(*pt3++)*(*pt4++);
	      }
	    }
	  }

	} else if (kj < nkernel) {

	  for (idx = 0; idx<nstamps; idx++) {

	    j1 = max(0,(int)stamp_ypos[idx]-stamp_half_width);
	    j2 = min(ny,(int)stamp_ypos[idx]+stamp_half_width);
	    i1 = max(0,(int)stamp_xpos[idx]-stamp_half_width);
	    i2 = min(nx,(int)stamp_xpos[idx]+stamp_half_width);

	    for (j=j1; j<j2; j++) {
	      pBi1 = &tex0[i1+nx*j];
	      pBj1 = &tex0[i1+nx*j];
	      if (ext_basis[kj]) {
		pBj2 = &tex1[i1+c+nx*(j+d)];
	      } else {
		pBj2 = &tex0[i1+c+nx*(j+d)];
	      }
	      pt3 = &tex3[i1+nx*j];
	      pt4 = &tex4[i1+nx*j];
	      ppx = &px[i1];
	      pyj = py[j];
	      for (i=i1; i<i2; i++) {
		temp += (*ppx++)*pyj*(*pBi1++)*((*pBj2++)-(*pBj1++))*(*pt3++)*(*pt4++);
	      }
	    }
	  }

	} else {

	  for (idx = 0; idx<nstamps; idx++) {

	    j1 = max(0,(int)stamp_ypos[idx]-stamp_half_width);
	    j2 = min(ny,(int)stamp_ypos[idx]+stamp_half_width);
	    i1 = max(0,(int)stamp_xpos[idx]-stamp_half_width);
	    i2 = min(nx,(int)stamp_xpos[idx]+stamp_half_width);

	    for (j=j1; j<j2; j++) {
	      pBi1 = &tex0[i1+nx*j];
	      pt3 = &tex3[i1+nx*j];
	      pt4 = &tex4[i1+nx*j];
	      ppx = &px[i1];
	      pyj = py[j];
	      for (i=i1; i<i2; i++) {
		temp += (*ppx++)*pyj*(*pBi1++)*(*pt3++)*(*pt4++);
	      }
	    }
	  }

	} 


      } else if (ki < nkernel) {

	if (kj == 0) {

	  for (idx = 0; idx<nstamps; idx++) {

	    j1 = max(0,(int)stamp_ypos[idx]-stamp_half_width);
	    j2 = min(ny,(int)stamp_ypos[idx]+stamp_half_width);
	    i1 = max(0,(int)stamp_xpos[idx]-stamp_half_width);
	    i2 = min(nx,(int)stamp_xpos[idx]+stamp_half_width);

	    for (j=j1; j<j2; j++) {
	      pBi1 = &tex0[i1+nx*j];
	      if (ext_basis[ki]) {
		pBi2 = &tex1[i1+a+nx*(j+b)];
	      } else {
		pBi2 = &tex0[i1+a+nx*(j+b)];
	      }
	      pBj1 = &tex0[i1+nx*j];
	      pt3 = &tex3[i1+nx*j];
	      pt4 = &tex4[i1+nx*j];
	      ppx = &px[i1];
	      pyj = py[j];
	      for (i=i1; i<i2; i++) {
		temp += (*ppx++)*pyj*((*pBi2++)-(*pBi1++))*(*pBj1++)*(*pt3++)*(*pt4++);
	      }
	    }
	  }

	} else if (kj < nkernel) {

	  for (idx = 0; idx<nstamps; idx++) {

	    j1 = max(0,(int)stamp_ypos[idx]-stamp_half_width);
	    j2 = min(ny,(int)stamp_ypos[idx]+stamp_half_width);
	    i1 = max(0,(int)stamp_xpos[idx]-stamp_half_width);
	    i2 = min(nx,(int)stamp_xpos[idx]+stamp_half_width);

	    for (j=j1; j<j2; j++) {
	      pBi1 = &tex0[i1+nx*j];
	      pBj1 = &tex0[i1+nx*j];
	      if (ext_basis[ki]) {
		pBi2 = &tex1[i1+a+nx*(j+b)];
	      } else {
		pBi2 = &tex0[i1+a+nx*(j+b)];
	      }
	      if (ext_basis[kj]) {
		pBj2 = &tex1[i1+c+nx*(j+d)];
	      } else {
		pBj2 = &tex0[i1+c+nx*(j+d)];
	      }
	      pt3 = &tex3[i1+nx*j];
	      pt4 = &tex4[i1+nx*j];
	      ppx = &px[i1];
	      pyj = py[j];
	      for (i=i1; i<i2; i++) {
		temp += (*ppx++)*pyj*((*pBi2++)-(*pBi1++))*((*pBj2++)-(*pBj1++))*(*pt3++)*(*pt4++);
	      }
	    }
	  }

	} else {

	  for (idx = 0; idx<nstamps; idx++) {

	    j1 = max(0,(int)stamp_ypos[idx]-stamp_half_width);
	    j2 = min(ny,(int)stamp_ypos[idx]+stamp_half_width);
	    i1 = max(0,(int)stamp_xpos[idx]-stamp_half_width);
	    i2 = min(nx,(int)stamp_xpos[idx]+stamp_half_width);

	    for (j=j1; j<j2; j++) {
	      pBi1 = &tex0[i1+nx*j];
	      if (ext_basis[ki]) {
		pBi2 = &tex1[i1+a+nx*(j+b)];
	      } else {
		pBi2 = &tex0[i1+a+nx*(j+b)];
	      }
	      pt3 = &tex3[i1+nx*j];
	      pt4 = &tex4[i1+nx*j];
	      ppx = &px[i1];
	      pyj = py[j];
	      for (i=i1; i<i2; i++) {
		temp += (*ppx++)*pyj*((*pBi2++)-(*pBi1++))*(*pt3++)*(*pt4++);
	      }
	    }
	  }

	}


      } else {

	if (kj == 0) {

	  for (idx = 0; idx<nstamps; idx++) {

	    j1 = max(0,(int)stamp_ypos[idx]-stamp_half_width);
	    j2 = min(ny,(int)stamp_ypos[idx]+stamp_half_width);
	    i1 = max(0,(int)stamp_xpos[idx]-stamp_half_width);
	    i2 = min(nx,(int)stamp_xpos[idx]+stamp_half_width);

	    for (j=j1; j<j2; j++) {
	      pBj1 = &tex0[i1+nx*j];
	      pt3 = &tex3[i1+nx*j];
	      pt4 = &tex4[i1+nx*j];
	      ppx = &px[i1];
	      pyj = py[j];
	      for (i=i1; i<i2; i++) {
		temp += (*ppx++)*pyj*(*pBj1++)*(*pt3++)*(*pt4++);
	      }
	    }
	  }

	} else if (kj < nkernel) {

	  for (idx = 0; idx<nstamps; idx++) {

	    j1 = max(0,(int)stamp_ypos[idx]-stamp_half_width);
	    j2 = min(ny,(int)stamp_ypos[idx]+stamp_half_width);
	    i1 = max(0,(int)stamp_xpos[idx]-stamp_half_width);
	    i2 = min(nx,(int)stamp_xpos[idx]+stamp_half_width);

	    for (j=j1; j<j2; j++) {
	      pBj1 = &tex0[i1+nx*j];
	      if (ext_basis[kj]) {
		pBj2 = &tex1[i1+c+nx*(j+d)];
	      } else {
		pBj2 = &tex0[i1+c+nx*(j+d)];
	      }
	      pt3 = &tex3[i1+nx*j];
	      pt4 = &tex4[i1+nx*j];
	      ppx = &px[i1];
	      pyj = py[j];
	      for (i=i1; i<i2; i++) {
		temp += (*ppx++)*pyj*((*pBj2++)-(*pBj1++))*(*pt3++)*(*pt4++);
	      }
	    }
	  }

	} else {

	  for (idx = 0; idx<nstamps; idx++) {

	    j1 = max(0,(int)stamp_ypos[idx]-stamp_half_width);
	    j2 = min(ny,(int)stamp_ypos[idx]+stamp_half_width);
	    i1 = max(0,(int)stamp_xpos[idx]-stamp_half_width);
	    i2 = min(nx,(int)stamp_xpos[idx]+stamp_half_width);

	    for (j=j1; j<j2; j++) {
	      pt3 = &tex3[i1+nx*j];
	      pt4 = &tex4[i1+nx*j];
	      ppx = &px[i1];
	      pyj = py[j];
	      for (i=i1; i<i2; i++) {
		temp += (*ppx++)*pyj*(*pt3++)*(*pt4++);
	      }
	    }
	  }
	}

      }


      H[blockIdx+gridDimx*blockIdy] = temp;
      H[blockIdy+gridDimx*blockIdx] = temp;
      if ((ki>0) && (ki<nkernel) && (kj>0) && (kj<nkernel)) {
	idx1 = np + (ki-1)*ns;
	idy1 = np + (kj-1)*ns;
	H[(idx1+idy0)+gridDimx*(idy1+idx0)] = temp;
	H[(idy1+idx0)+gridDimx*(idx1+idy0)] = temp;
      }

    }
  }

}
