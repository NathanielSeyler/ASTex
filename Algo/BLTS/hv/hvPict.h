/*******************************************************************************
* ASTex:                                                                       *
* Copyright (C) IGG Group, ICube, University of Strasbourg, France             *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: https://astex-icube.github.io                                      *
* Contact information: astex@icube.unistra.fr                                  *
*                                                                              *
*******************************************************************************/



// hvPicture.h: interface for the picture class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PICT_H__098453F0_1C38_49E9_A6F4_AABF90AA55E8__INCLUDED_)
#define AFX_PICT_H__098453F0_1C38_49E9_A6F4_AABF90AA55E8__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include "hvField2.h"
#include "hvColor.h"
#include "hvList.h"
#include "hvVec2.h"

namespace hview {

template <class T> class hvPictRGB ;
template <class T> class hvPictRGBA ;
class hvBitmap;

enum hvPictComponent { HV_RED, HV_GREEN, HV_BLUE, HV_ALPHA, HV_RGB, HV_LUMINANCE } ;
enum hvPictMask { HV_GAUSS, HV_EDGE_VERT, HV_EDGE_HORIZ, HV_EDGE_DIAG1, HV_EDGE_DIAG2, HV_DESCR_HORIZ, HV_DESCR_VERT, HV_DESCR_DIAG1, HV_DESCR_DIAG2, HV_DESCR_DOTS } ;

////////////////////////////////////////////////////////////
template <class T> class hvPict : public hvField2< T >  // T is a scalar type
////////////////////////////////////////////////////////////
{
public:
	hvPict<T>() : hvField2< T >(),hvArray2< T >() { }
	hvPict<T>(int sx, int sy, T nil) : hvField2< T >(sx, sy, nil),hvArray2< T >(sx, sy, nil) { }
	template <class X> hvPict<T>(const hvArray2<X> &pict)
	{
		hvField2< T >::reset(pict.sizeX(), pict.sizeY(), T(0));
		int i, j;
		for (i = 0; i<pict.sizeX(); i++) for (j = 0; j<pict.sizeY(); j++)
		{
			update(i, j, T(pict.get(i, j)));
		}
	}
	template <class X> void clone(const hvPict<X> &pict)
	{
		hvField2< T >::reset(pict.sizeX(), pict.sizeY(), T(0));
		int i,j;
		for (i=0; i<pict.sizeX(); i++) for (j=0; j<pict.sizeY(); j++)
		{
			update(i,j,T(pict.get(i,j)));
		}
	}
	template <class X> void clone(const hvPict<X> &pict, int x, int y, int sx, int sy)
	{
		hvField2< T >::reset(sx-x+1, sy-y+1, T(0));
		int i,j;
		for (i=x; i<=sx; i++) for (j=y; j<=sy; j++)
		{
			update(i-x,j-y,T(pict.get(i,j)));
		}
	}
	hvPict<T>(const hvBitmap &pict, T va, T vb, int x, int y, int sx, int sy): hvField2< T >(sx-x+1, sy-y+1, T(0)), hvArray2< T >(sx-x+1, sy-y+1, T(0))
	{
		int i,j;
		for (i=x; i<=sx; i++) for (j=y; j<=sy; j++)
		{
			if (pict.get(i,j)) update(i-x,j-y,va); else update(i-x,j-y,vb);
		}
	}
	hvPict<T>(const hvBitmap &pict, int border, T scal): hvField2< T >(pict.sizeX(), pict.sizeY(), T(0)), hvArray2< T >(pict.sizeX(), pict.sizeY(), T(0))
	{
		hvPict<int> pi(pict.sizeX(), pict.sizeY(), 0);
		hvBitmap bm; bm=pict;
		int i,j;
		int count =1;
		bool cont=true;
		do {
			hvBitmap bb; bb.erosion(bm,3,3); 
			hvBitmap back; back=bb;
			//printf("thinning %d, %d\n", count, bm.count());
			~bb; bb &= bm;
			for (i=0; i<pict.sizeX(); i++) for (j=0; j<pict.sizeY(); j++)
			{
				if (bb.get(i,j)) pi.update(i,j,count);
			}
			count++;
			if (count>border) count=border;
			bm = back;
			//if (!cont) { for (i=0; i<pict.sizeX(); i++) for (j=0; j<pict.sizeY(); j++) if (bb.get(i,j)) pi.update(i,j,count); }
		} while (bm.count()>0);
		for (i=0; i<pict.sizeX(); i++) for (j=0; j<pict.sizeY(); j++)
		{
			this->update(i,j,T((double)pi.get(i,j)/(double)border*(double)scal));
		}

	}
	
	hvPict<T>(bool loop, const hvBitmap &pict, int border, T scal): hvField2< T >(pict.sizeX(), pict.sizeY(), T(0)), hvArray2< T >(pict.sizeX(), pict.sizeY(), T(0))
	{
		hvPict<int> pi(pict.sizeX(), pict.sizeY(), border);
		hvBitmap bm; bm=pict;
		int i,j;
		int count =1;
		int upd;
		hvVec2<int> bmin, bmax;
		pict.box(bmin, bmax);
		int minx=bmin.X(), miny=bmin.Y(), maxx=bmax.X(), maxy=bmax.Y();
		bool cont=true;
		do {
			hvBitmap bb; bb.erosionTorus(bm,3,3, minx,miny, maxx, maxy); 
			hvBitmap back; back=bb;
			//printf("thinning %d, %d\n", count, bm.count());
			~bb; bb &= bm;
			upd=0;
			for (i=minx; i<=maxx; i++) for (j=miny; j<=maxy; j++)
			{
				if (bb.get(i,j)) { upd++; pi.update(i,j,count); }
			}
			count++;
			if (count>border) cont=false;
			bm = back;
			//if (!cont) { for (i=0; i<pict.sizeX(); i++) for (j=0; j<pict.sizeY(); j++) if (bb.get(i,j)) pi.update(i,j,count); }
		} while (cont && upd>0);
		for (i=minx; i<=maxx; i++) for (j=miny; j<=maxy; j++)
		{
			this->update(i,j,T((double)pi.get(i,j)/(double)border*(double)scal));
		}

	}

	hvPict<T>(const hvBitmap &pict, int border, T scal, int kernel, int startx, int starty, int endx, int endy): hvField2< T >(pict.sizeX(), pict.sizeY(), T(0)), hvArray2< T >(pict.sizeX(), pict.sizeY(), T(0))
	{
		this->reset(pict.sizeX(), pict.sizeY(), T(0));
		hvPict<int> pi(pict.sizeX(), pict.sizeY(), 0);
		hvBitmap bm; bm=pict;
		int i,j;
		int count =1;
		bool cont=true;
		do {
			hvBitmap bb; bb.erosion(bm,kernel,kernel, startx, starty,endx, endy); 
			hvBitmap back; back=bb;
			//printf("thinning %d\n", count);
			~bb; bb &= bm;
			for (i=startx; i<=endx; i++) for (j=starty; j<=endy; j++)
			{
				if (i>=0 && i<bb.sizeX() && j>=0 && j<bb.sizeY()) if (bb.get(i,j)) pi.update(i,j,count);
			}
			count++;
			if (count>border) count=border;
			bm = back;
			//if (!cont) { for (i=0; i<pict.sizeX(); i++) for (j=0; j<pict.sizeY(); j++) if (bb.get(i,j)) pi.update(i,j,count); }
		} while (bm.count()>0);
		for (i=0; i<pict.sizeX(); i++) for (j=0; j<pict.sizeY(); j++)
		{
			this->update(i,j,T((double)pi.get(i,j)/(double)border*(double)scal));
		}

	}
	template <class X, class Y> hvPict<T>(const hvPict<X> &pict, Y scal, int x, int y, int sx, int sy): hvField2< T >(sx-x+1, sy-y+1, T(0)), hvArray2< T >(sx-x+1, sy-y+1, T(0))
	{
		this->reset(pict.sizeX(), pict.sizeY(), T(0));
		int i,j;
		for (i=0; i<=sx; i++) for (j=0; j<=sy; j++)
		{
			this->update(i,j,T(Y(pict.get(i,j))*scal));
		}
	}
	template <class X> hvPict<T>(const hvPict<X> &pict, double max, double loga, double scal, int x, int y, int sx, int sy) : hvField2< T >(sx - x + 1, sy - y + 1, T(0)), hvArray2< T >(sx - x + 1, sy - y + 1, T(0))
	{
		this->reset(pict.sizeX(), pict.sizeY(), T(0));
		int i, j;
		for (i = 0; i <= sx; i++) for (j = 0; j <= sy; j++)
		{
			X val = pict.get(i, j);
			double vv = log(1.0 + loga*(double)val / max) / log(loga + 1.0);
			if (vv>1.0) vv = 1.0;
			this->update(i, j, T(vv*scal));
		}
	}

	template <class X, class Y> hvPict<T>(const hvPictRGB<X> &pict, Y scal, int x, int y, int sx, int sy);
	template <class X, class Y> void convert(const hvPictRGB<X> &pict, hvPictComponent cc, Y scal);
	template <class X, class Y> hvPict<T>(const hvPictRGB<X> &pict, hvPictComponent cc, Y scal, int x, int y, int sx, int sy);
	template <class X, class Y> hvPict<T>(const hvPictRGBA<X> &pict, Y scal, int x, int y, int sx, int sy);
	template <class X, class Y> hvPict<T>(const hvPictRGBA<X> &pict, hvPictComponent cc, Y scal, int x, int y, int sx, int sy);
	
	template <class X, class Y> hvPict<T>(hvArray1<X> *ft, int size, X loga, Y scal): hvField2< T >(size, size, T(0)), hvArray2< T >(size, size, T(0))
	{
		int i,j;
		X max = X(0);
		for (i=0; i<size; i++) for (j=0; j<size; j++)
		{
			if (i!=size/2 || j!=size/2)
			{
				X val = ft->get(j+i*size);
				if (max<val) max=val;
			}
		}
		for (i=0; i<size; i++) for (j=0; j<size; j++)
		{
			X val = ft->get(j+i*size);
			val = (X)log(1.0+loga*val/max)/(X)log(loga+1.0);
			if (val>X(1)) val=X(1);
			update(j,i, T(Y(val)*scal));
		}
	}
	template <class X, class Y> hvPict<T>(hvArray1<hvPair<X, X> > *ft, int size, X loga, Y scal, bool amplitude=false): hvField2< T >(size, size, T(0)), hvArray2< T >(size, size, T(0))
	{
		int i,j;
		X max = X(0);
		for (i=0; i<size; i++) for (j=0; j<size; j++)
		{
			if (i!=size/2 || j!=size/2)
			{
				hvPair<X, X> vv = ft->get(j+i*size);
				X val;
				if (amplitude) val = sqrt(vv.getLeft()*vv.getLeft()+vv.getRight()*vv.getRight());
				else val = vv.getLeft()*vv.getLeft()+vv.getRight()*vv.getRight();
				if (max<val) max=val;
			}
		}
		for (i=0; i<size; i++) for (j=0; j<size; j++)
		{
			hvPair<X, X> vv = ft->get(j+i*size);
			X val;
			if (amplitude) val = sqrt(vv.getLeft()*vv.getLeft()+vv.getRight()*vv.getRight());
			else val = vv.getLeft()*vv.getLeft()+vv.getRight()*vv.getRight();
			val = (X)log(1.0+loga*val/max)/(X)log(loga+1.0);
			update(j,i, T(Y(val)*scal));
		}
	}
	// devides size by factor of 2^(nrec+1), nrec=0 means factor 2, nrec=1 means factor 4, etc.
	void shrink(hvPict<T> *p, int nrec=0)
	{
		if (nrec>0) { hvPict<T> source; source.hvField2<T>::shrink(p); this->shrink(&source, nrec-1); }
		else hvField2<T>::shrink(p);
	}
	// multiplies size by factor of 2^(nrec+1), nrec=0 means factor 2, nrec=1 means factor 4, etc.
	void enlarge(hvPict<T> *p, int nrec=0)
	{
		if (nrec>0) { hvPict<T> source; source.hvField2<T>::enlarge(p); this->enlarge(&source, nrec-1); }
		else hvField2<T>::enlarge(p);
	}
	// applies a deblurring mask of size sx,sy that must be odd numbers
	void deblur(int niter, int sx, int sy, double scal, double min, double max)
	{
		hvField2<T> source; source.clone(this);
		hvField2<T>:: template deblurs<double>(&source, niter, sx, sy, scal, min, max);
	}
	// applies a blurring mask of size sx,sy
	// scal is the normalization factor, usually sx*sy
	void blur(T scal, int sx, int sy)
	{
		hvField2<T> source; source.clone(this);
		hvField2<T>:: template blur<double>(&source, sx, sy,(double)(1.0/(double)scal));
	}
	void gaussianBlur(int sx, int sy)
	{
		hvField2<T> source; source.clone(this);
		hvField2<T>::template gaussianBlur<double>(&source, sx, sy);
	}
	void gaussianBlur(int sx, int sy, double sig)
	{
		hvField2<T> source; source.clone(this);
		hvField2<T>:: template gaussianBlur<double>(&source, sx, sy, sig);
	}
	void gaborFilter(int sx, int sy, double theta, double freq, double sigma, double phase)
	{
		hvField2<T> source; source.clone(this);
		hvField2<T>:: template gaborFilter<double>(&source, sx, sy, theta, freq, sigma, phase);
	}
	void anisoGaborFilter(int sx, int sy, double alpha, double theta, double freq, double sigmax, double sigmay, double phase)
	{
		hvField2<T> source; source.clone(this);
		hvField2<T>:: template anisoGaborFilter<double>(&source, sx, sy, alpha, theta, freq, sigmax, sigmay, phase);
	}
	// alpha1,2 must be between -PI/4 and 3PI/4 
	void bandpassFilter(hvPict<T> &res, int pow2, double f1, double f2, double alpha1, double alpha2, double offset, double scale)
	{
		int i,j, ii,jj;
		int NN=1<<pow2;
		res.hvField2< T >::reset(this->sizeX(),this->sizeY(),T(0));
		for (i=0; i<this->sizeX(); i++)
		{
			//if (i%10==0) { printf("."); fflush(stdout); }
			for (j=0; j<this->sizeY(); j++)
			{
				hvField2<double> sub(NN,NN,0.0);
				int xx,yy;
				for (ii=0; ii<NN; ii++) for (jj=0; jj<NN;jj++) 
				{
					xx = i+ii-NN/2; if (xx<0) xx+= this->sizeX(); else if (xx>= this->sizeX()) xx-= this->sizeX();
					yy = j+jj-NN/2; if (yy<0) yy+= this->sizeY(); else if (yy>= this->sizeY()) yy-= this->sizeY();
					sub.update(ii,jj,(double)this->get(xx,yy));
				}
				hvArray1<hvPair<double,double> > *ft = sub.fft(true, scale, offset, pow2);
				for (ii=0; ii<NN; ii++) for (jj=0; jj<NN; jj++)
				{
							double c = ft->get(ii+jj*NN).mod();
							double phase = ft->get(ii+jj*NN).phase();
							double fx = (double)(ii-NN/2)/(double)NN;
							double fy = (double)(jj-NN/2)/(double)NN;
							double ff = sqrt(fx*fx+fy*fy);
							if (ff<f1 || ff>f2 ) c=0.0;
							double alpha=hvPair<double,double>(fx,fy).phase();
							if (alpha<-M_PI/4.0) alpha+=M_PI;
							else if (alpha>3.0*M_PI/4.0) alpha-=M_PI;
							if (alpha<alpha1 || alpha>alpha2) c=0.0;
							ft->update(ii+jj*NN, hvPair<double,double>( c*cos(phase), -c*sin(phase) ));
				}
				for (ii=0; ii<NN; ii++) hvArray1<double>::fft(*ft, pow2,1,ii*NN,false);
				for (ii=0; ii<NN; ii++) hvArray1<double>::fft(*ft, pow2,NN,ii,false);
				hvPair<double,double> c = ft->get(NN/2+NN/2*NN);
				double val=c.mod()*cos(c.phase());
				//printf("%d,%d= %g versus %g\n",i,j, (double)this->get(i,j), val);
				res.update(i,j,val*scale+offset);
				delete ft;
			}
		}
	}
	// alpha1,2 must be between -PI/4 and 3PI/4 
	void bandpassFilterBank(int n, hvPict<T> res[], int pow2, double f1[], double f2[], double alpha1[], double alpha2[], bool inout[], double offset, double scale)
	{
		int i,j, k, ii,jj;
		int NN=1<<pow2;
		for (i=0; i<n; i++) res[i].hvField2< T >::reset(this->sizeX(), this->sizeY(),T(0));
		for (i=0; i<this->sizeX(); i++)
		{
			//if (i%10==0) { printf("."); fflush(stdout); }
			for (j=0; j<this->sizeY(); j++)
			{
				hvField2<double> sub(NN,NN,0.0);
				int xx,yy;
				for (ii=0; ii<NN; ii++) for (jj=0; jj<NN;jj++) 
				{
					xx = i+ii-NN/2; if (xx<0) xx+= this->sizeX(); else if (xx>= this->sizeX()) xx-= this->sizeX();
					yy = j+jj-NN/2; if (yy<0) yy+= this->sizeY(); else if (yy>= this->sizeY()) yy-= this->sizeY();
					sub.update(ii,jj,(double)this->get(xx,yy));
				}
				hvArray1<hvPair<double,double> > *ft = sub.fft(true, scale, offset, pow2);
				hvArray1<hvPair<double,double> > *rft = new hvArray1<hvPair<double,double> >(NN*NN,hvPair<double,double>(0.0,0.0));
				for (k=0; k<n; k++)
				{
					for (ii=0; ii<NN; ii++) for (jj=0; jj<NN; jj++)
					{
								double c = ft->get(ii+jj*NN).mod();
								double phase = ft->get(ii+jj*NN).phase();
								double fx = (double)(ii-NN/2)/(double)NN;
								double fy = (double)(jj-NN/2)/(double)NN;
								double ff = sqrt(fx*fx+fy*fy);
								if (ff<f1[k] || ff>f2[k] ) c=0.0;
								double alpha=hvPair<double,double>(fx,fy).phase();
								if (alpha<-M_PI/4.0) alpha+=M_PI;
								else if (alpha>3.0*M_PI/4.0) alpha-=M_PI;
								if (inout[k]) { if (alpha<alpha1[k] || alpha>alpha2[k]) c=0.0; }
								else { if (alpha>=alpha1[k] && alpha<=alpha2[k]) c=0.0; }
								rft->update(ii+jj*NN, hvPair<double,double>( c*cos(phase), -c*sin(phase) ));
					}
					for (ii=0; ii<NN; ii++) hvArray1<double>::fft(*rft, pow2,1,ii*NN,false);
					for (ii=0; ii<NN; ii++) hvArray1<double>::fft(*rft, pow2,NN,ii,false);
					hvPair<double,double> c = rft->get(NN/2+NN/2*NN);
					double val=c.mod()*cos(c.phase());
					//printf("%d,%d= %g versus %g\n",i,j, (double)this->get(i,j), val);
					res[k].update(i,j,val*scale+offset);
				}
				delete ft;
				delete rft;
			}
		}
	}

	// bilateral filter
	void bilateral(double sigma, T scal, int sx, int sy)
	{
		hvPict<T> pp(this->sizeX(), this->sizeY(), scal);
		pp.clone(*this);
		this->bilateral(pp, sigma, scal, sx, sy);
	}
	void bilateral(const hvPict<T> &pia, double sigma, T scal, int sx, int sy)
	{
		hvField2< T >::reset(pia.sizeX(), pia.sizeY(), T(0));
		int i, j, ii, jj;
		for (i = 0; i<pia.sizeX(); i++) for (j = 0; j<pia.sizeY(); j++)
		{
			T cola = pia.get(i, j);
			double ra = (double)cola / (double)scal;
			double ww = 0.0;
			double resa = 0.0;
			for (ii = -sx / 2; ii <= sx / 2; ii++) for (jj = -sy / 2; jj <= sy / 2; jj++)
			{
				if (i + ii >= 0 && i + ii<pia.sizeX() && j + jj >= 0 && j + jj<pia.sizeY())
				{
					cola = pia.get(i + ii, j + jj);
					double va = (double)cola / (double)scal;
					double dist = (va - ra)*(va - ra);
					dist = exp(-dist*sigma*sigma);
					ww += dist;
					resa += dist*va;
				}
			}
			resa /= ww;
			update(i, j, T(resa*scal));
		}
	}

	void median(int ss)
	{
		hvField2<T> source; source.clone(this);
		hvField2<T>::median(&source, ss);
	}

	T maxValue() const
	{
		T m = this->get(0,0);
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			T v = this->get(i,j);
			if (v>m) m=v;
		}
		return m;
	}
	T maxValue(int x, int y, int dx, int dy) const
	{
		T m = this->get(0,0);
		int i,j;
		for (i=x;i<dx; i++) for (j=y; j<dy; j++)
		{
			T v = this->get(i,j);
			if (v>m) m=v;
		}
		return m;
	}
	T minValue() const
	{
		T m = this->get(0,0);
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			T v = this->get(i,j);
			if (v<m) m=v;
		}
		return m;
	}

	void normalize(T scal, double percent=1.0, double power=1.0) 
	{
		int i,j;
		double max=(double)maxValue()*percent;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			double v = (double)this->get(i,j)/max;
			if (v>1.0) v=1.0;
			this->update(i,j,T(pow(v,power)*(double)scal));
		}
	}
	void normalize(T scal, T maxval, double percent, double power) 
	{
		int i,j;
		double max=(double)maxval*percent;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			double v = (double)this->get(i,j)/max;
			if (v>1.0) v=1.0;
			this->update(i,j,T(pow(v,power)*(double)scal));
		}
	}

	void rescale(T nmin, T nmax) 
	{
		int i,j;
		T min, max;
		this->minmax(min, max);
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			double v = (double)nmin+((double)nmax-(double)nmin)*((double)this->get(i,j)-(double)min)/((double)max-(double)min);
			this->update(i,j,T(v));
		}
	}
	void reverse(T norm) 
	{
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			this->update(i,j,norm-this->get(i,j));
		}
	}

	void clamp() 
	{
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			T v = this->get(i,j);
			if (v<T(0)) v=T(0);
			this->update(i,j,v);
		}
	}

	void gamma(T scal, double power) 
	{
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			T v = this->get(i,j);
			v = (T)((double)scal*pow((double)v/(double)scal,power));
			this->update(i,j,v);
		}
	}
	void contrast(T scal, double power) 
	{
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			double v = (double)this->get(i,j);
			v /= (double)scal;
			if (v>=0.5) v = pow((v-0.5)*2.0,1.0/power)/2.0+0.5;
			else v = pow(v*2.0, power)/2.0;
			this->update(i,j,T(v*(double)scal));
		}
	}
	double avg() const 
	{
		int i,j;
		double v=0.0;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			v += (double)this->get(i,j);
		}
		return v/(double)(this->sizeX()*this->sizeY());
	}
	double avg(const hvBitmap &mask) const 
	{
		int i,j;
		double v=0.0;
		int n=0;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			if (mask.get(i,j)) { n++; v += (double)this->get(i,j); }
		}
		return v/(double)n;
	}
	double avgrange(const hvBitmap &mask, T min, T max) const 
	{
		int i,j;
		double v=0.0;
		int n=0;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			T val = this->get(i,j);
			if (mask.get(i,j) && val>=min && val<=max) { n++; v += (double)val; }
		}
		return v/(double)n;
	}
	double avgstat(T &min, T &max, double &var) const 
	{
		int i,j;
		double v=0.0, v2=0.0;
		min=this->get(0,0); max=this->get(0,0);
		int n=this->sizeX()*this->sizeY();
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			T val = this->get(i,j);
			v += (double)val;
			v2 += (double)val*(double)val;
			if (val<min) min=val;
			if (val>max) max=val;
		}
		v /= (double)n;
		var = v2/(double)n-v*v;
		return v;
	}

	template <class X> void histogramm(hvArray1<int> &tab, X min, X max) const
	{
		int i,j;
		for (i=0; i<tab.size(); i++) tab.update(i,0);
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			T v = this->get(i,j);
			X ind = (X(v)-min)/(max-min);
			int pos = (int)(ind*X(tab.size()-1));
			tab.update(pos,tab.get(pos)+1);
		}
	}
	template <class U> void blendRect(int px, int py, int x, int y, int sx, int sy, const hvPict<T> &example, const hvPict<U> &alpha, U scal, double power, const hvBitmap &mask)
	{
		int i,j;
		for (i=0; i<sx; i++) for (j=0; j<sy; j++)
		{
			if (mask.get(x+i,y+j))
			{
				if (px+i>=0 && px+i<this->sizeX() && py+j>=0 && py+j<this->sizeY())
				{
					double coeff = pow((double)alpha.get(x+i,y+j)/(double)scal, power);
					T col = this->get(px+i,py+j);
					T colex = example.get(x+i,y+j);
					T colres = T((double)colex*coeff+(double)col*(1.0-coeff) );
					this->update(px+i,py+j,colres);
				}
			}
		}
	}
	/*
	template <class X> void gnoise(X scal, X offset, double fx, double fy)
	{
		int i,j;
		for (i=0;i<sizeX(); i++) for (j=0; j<sizeY(); j++)
		{
			X v = (X)hvNoise::gnoise((double)i*(double)fx+23.987, (double)j*(double)fy-6.3327, 42.12345)*scal+offset;
			update(i,j,T(v));
		}
	}
	template <class X> void snoise(X scal, X offset, double fx, double fy)
	{
		int i,j;
		for (i=0;i<sizeX(); i++) for (j=0; j<sizeY(); j++)
		{
			X v = (X)hvNoise::snoise2((double)i*(double)fx+23.987, (double)j*(double)fy-6.3327)*scal+offset;
			update(i,j,T(v));
		}
	}
	*/
	bool localMaximum(int x, int y, int sx, int sy) const
	{
		int i,j;
		for(i=x-sx; i<=x+sx; i++) for (j=y-sy; j<=y+sy; j++)
		{
			if (i>=0 && i<this->sizeX() && j>=0 && j<this->sizeY() && (x!=i || y!=j))
			{
				if (this->get(i,j)>=this->get(x,y)) return false;
			}
		}
		return true;
	}
	void minmax(const hvBitmap &mask, T &min, T &max) const
	{
		int i,j;
		bool first=true;
		for(i=0; i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			if (mask.get(i,j))
			{
				if (first) { min=this->get(i,j); max=this->get(i,j); first=false; }
				else
				{
					T v = this->get(i,j);
					if (v<min) min=v;
					if (v>max) max=v;
					//printf("%d,%d=%d\n", i,j,v);
				}
			}
		}
	}
	void minmax(T &min, T &max) const
	{
		int i,j;
		bool first=true;
		for(i=0; i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
				if (first) { min=this->get(i,j); max=this->get(i,j); first=false; }
				else
				{
					T v = this->get(i,j);
					if (v<min) min=v;
					if (v>max) max=v;
					//printf("%d,%d=%d\n", i,j,v);
				}
		}
	}

	float mse(const hvPict<T> &pp, float scal)
	{
		int i,j;
		float res=0;
		for (i=0; i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			float v1 = (float)this->get(i,j)/scal; 
			float v2 = (float)pp.get(i,j)/scal;
			res += (v1-v2)*(v1-v2);
		}
		return res;
	}

	void squaredDifference(int px, int py, int dx, int dy, const hvPictRGB<unsigned char> &pia, int ix, int iy, const hvPictRGB<unsigned char> &pib);

	double minPathH(int ystart, int yend, int mm, int ypath[], double *herr=0)
	{
		int i,j;
		hvArray2<double> sumdiff(this->sizeX(),this->sizeY(),0.0);
		for (j=0; j<this->sizeY(); j++) sumdiff.update(0,j,(double)this->get(0,j));
		hvArray2<int> prev(this->sizeX(),this->sizeY(),0);
		// compute min cut path
		for (i=1; i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			int mink;
			double sumdiffmin=0.0;
			bool first=true;
			for (int k=-mm; k<=mm; k++)
			{
				if (j+k>=0 && j+k<this->sizeY())
					if (sumdiff.get(i-1,j+k)<sumdiffmin || first) { sumdiffmin=sumdiff.get(i-1,j+k); mink=k; first=false; }
			}
			sumdiff.update(i,j, sumdiff.get(i-1,j+mink)+(double)this->get(i,j));
			prev.update(i,j,mink);
			/*
			if (j==0)
			{
				if (sumdiff.get(i-1,j)<sumdiff.get(i-1,j+1)) 
				{
					sumdiff.update(i,j, sumdiff.get(i-1,j)+(double)this->get(i,j));
					prev.update(i,j,0);
					//std::cout<<i<<","<<j<<"="<<sumdiff.get(i,j)<<";"<<sqdiff.get(i,j)<<"\n";
				}
				else
				{
					sumdiff.update(i,j, sumdiff.get(i-1,j+1)+(double)this->get(i,j));
					prev.update(i,j,1);
					//std::cout<<i<<","<<j<<"="<<sumdiff.get(i,j)<<";"<<sqdiff.get(i,j)<<"\n";
				}
			}
			else if (j==this->sizeY()-1)
			{
				if (sumdiff.get(i-1,j)<sumdiff.get(i-1,j-1)) 
				{
					sumdiff.update(i,j, sumdiff.get(i-1,j)+(double)this->get(i,j));
					prev.update(i,j,0);
					//std::cout<<i<<","<<j<<"="<<sumdiff.get(i,j)<<";"<<sqdiff.get(i,j)<<"\n";
				}
				else
				{
					sumdiff.update(i,j, sumdiff.get(i-1,j-1)+(double)this->get(i,j));
					prev.update(i,j,-1);
					//std::cout<<i<<","<<j<<"="<<sumdiff.get(i,j)<<";"<<sqdiff.get(i,j)<<"\n";
				}
			}
			else
			{
				if (sumdiff.get(i-1,j)<=sumdiff.get(i-1,j-1) && sumdiff.get(i-1,j)<=sumdiff.get(i-1,j+1)) 
				{
					sumdiff.update(i,j, sumdiff.get(i-1,j)+(double)this->get(i,j));
					prev.update(i,j,0);
					//std::cout<<i<<","<<j<<"="<<sumdiff.get(i,j)<<";"<<sqdiff.get(i,j)<<"\n";
				}
				else if (sumdiff.get(i-1,j-1)<=sumdiff.get(i-1,j) && sumdiff.get(i-1,j-1)<=sumdiff.get(i-1,j+1))
				{
					sumdiff.update(i,j, sumdiff.get(i-1,j-1)+(double)this->get(i,j));
					prev.update(i,j,-1);
					//std::cout<<i<<","<<j<<"="<<sumdiff.get(i,j)<<";"<<sqdiff.get(i,j)<<"\n";
				}
				else
				{
					sumdiff.update(i,j, sumdiff.get(i-1,j+1)+(double)this->get(i,j));
					prev.update(i,j,1);
					//std::cout<<i<<","<<j<<"="<<sumdiff.get(i,j)<<";"<<sqdiff.get(i,j)<<"\n";
				}
			}
			*/
		}
		if (ystart<0) ystart=0;
		if (yend>=this->sizeY()) yend=this->sizeY()-1;
		int miny=ystart; 
		double minv=sumdiff.get(this->sizeX()-1,ystart);
		for (j=ystart+1; j<=yend; j++) if (sumdiff.get(this->sizeX()-1,j)<minv) { minv=sumdiff.get(this->sizeX()-1,j); miny=j; }
		int posy=miny;
		double total = 0.0;
		if (herr != 0) { herr[0] = 0.0;  herr[1] = 0.0; }
		for (i=0; i<this->sizeX(); i++)
		{
			//std::cout<<posy<<"\n";
			double ee= this->get(this->sizeX() - 1 - i, posy);
			total += ee;
			if (herr != 0) { if (i < this->sizeX() / 2) herr[0] += ee; else herr[1] += ee; }
			ypath[this->sizeX()-1-i]=posy;
			posy += prev.get(this->sizeX()-1-i,posy);
		}
		if (herr != 0) { herr[0] /= (double)(this->sizeX()/2);  herr[1] = (double)(this->sizeX() / 2); }
		return total/(double)(this->sizeX());
	}
	double minPathV(int xstart, int xend, int mm, int xpath[], double *herr = 0)
	{
		int i,j;
		hvArray2<double> sumdiff(this->sizeX(),this->sizeY(),0.0);
		for (i=0; i<this->sizeX(); i++) sumdiff.update(i,0,(double)this->get(i,0));
		hvArray2<int> prev(this->sizeX(),this->sizeY(),0);
		// compute min cut path
		for (j=1; j<this->sizeY(); j++) for (i=0; i<this->sizeX(); i++)
		{
			int mink;
			double sumdiffmin=0;
			bool first=true;
			for (int k=-mm; k<=mm; k++)
			{
				if (i+k>=0 && i+k<this->sizeX())
					if (sumdiff.get(i+k,j-1)<sumdiffmin || first) { sumdiffmin=sumdiff.get(i+k,j-1); mink=k; first=false; }
			}
			sumdiff.update(i,j, sumdiff.get(i+mink,j-1)+(double)this->get(i,j));
			prev.update(i,j,mink);
			/*
			if (i==0)
			{
				if (sumdiff.get(i,j-1)<sumdiff.get(i+1,j-1)) 
				{
					sumdiff.update(i,j, sumdiff.get(i,j-1)+(double)this->get(i,j));
					prev.update(i,j,0);
					//std::cout<<i<<","<<j<<"="<<sumdiff.get(i,j)<<";"<<sqdiff.get(i,j)<<"\n";
				}
				else
				{
					sumdiff.update(i,j, sumdiff.get(i+1,j-1)+(double)this->get(i,j));
					prev.update(i,j,1);
					//std::cout<<i<<","<<j<<"="<<sumdiff.get(i,j)<<";"<<sqdiff.get(i,j)<<"\n";
				}
			}
			else if (i==this->sizeX()-1)
			{
				if (sumdiff.get(i,j-1)<sumdiff.get(i-1,j-1)) 
				{
					sumdiff.update(i,j, sumdiff.get(i,j-1)+(double)this->get(i,j));
					prev.update(i,j,0);
					//std::cout<<i<<","<<j<<"="<<sumdiff.get(i,j)<<";"<<sqdiff.get(i,j)<<"\n";
				}
				else
				{
					sumdiff.update(i,j, sumdiff.get(i-1,j-1)+(double)this->get(i,j));
					prev.update(i,j,-1);
					//std::cout<<i<<","<<j<<"="<<sumdiff.get(i,j)<<";"<<sqdiff.get(i,j)<<"\n";
				}
			}
			else
			{
				if (sumdiff.get(i,j-1)<=sumdiff.get(i-1,j-1) && sumdiff.get(i,j-1)<=sumdiff.get(i+1,j-1)) 
				{
					sumdiff.update(i,j, sumdiff.get(i,j-1)+(double)this->get(i,j));
					prev.update(i,j,0);
					//std::cout<<i<<","<<j<<"="<<sumdiff.get(i,j)<<";"<<sqdiff.get(i,j)<<"\n";
				}
				else if (sumdiff.get(i-1,j-1)<=sumdiff.get(i,j-1) && sumdiff.get(i-1,j-1)<=sumdiff.get(i+1,j-1))
				{
					sumdiff.update(i,j, sumdiff.get(i-1,j-1)+(double)this->get(i,j));
					prev.update(i,j,-1);
					//std::cout<<i<<","<<j<<"="<<sumdiff.get(i,j)<<";"<<sqdiff.get(i,j)<<"\n";
				}
				else
				{
					sumdiff.update(i,j, sumdiff.get(i+1,j-1)+(double)this->get(i,j));
					prev.update(i,j,1);
					//std::cout<<i<<","<<j<<"="<<sumdiff.get(i,j)<<";"<<sqdiff.get(i,j)<<"\n";
				}
			}
			*/
		}
		if (xstart<0) xstart=0;
		if (xend>=this->sizeX()) xend=this->sizeX()-1;
		int minx=xstart; 
		double minv=sumdiff.get(xstart,this->sizeY()-1);
		for (i=xstart+1; i<=xend; i++) if (sumdiff.get(i,this->sizeY()-1)<minv) { minv=sumdiff.get(i,this->sizeY()-1); minx=i; }
		int posx=minx;
		double total = 0.0;
		if (herr != 0) { herr[0] = 0.0;  herr[1] = 0.0; }
		for (j=0; j<this->sizeY(); j++)
		{
			//std::cout<<posx<<"\n";
			double ee = this->get(posx, this->sizeY() - 1 - j);
			total += ee;
			if (herr != 0) { if (j < this->sizeY() / 2) herr[0] += ee; else herr[1] += ee; }
			xpath[this->sizeY()-1-j]=posx;
			posx += prev.get(posx, this->sizeY()-1-j);
		}
		if (herr != 0) { herr[0] /= (double)(this->sizeY() / 2);  herr[1] = (double)(this->sizeY() / 2); }
		return total / (double)(this->sizeY());
	}
	void wft(const hvBitmap &mask, int size, float sc, const std::vector<hvVec2<int> > &ff, hvPict<float> *histo)
	{
		if (histo->sizeY()!=ff.size()) { hvFatal("histo sizeY does not match number of frequencies!"); }
		//T min,max;
		int nn=0;
		//minmax(mask, min,max);
		//printf("minmax=%d,%d\n", min,max);
		hvField2<float> win(size,size,0.0f);
		hvArray1<hvPair<float,float> > *fft;
		int i,j, ii, jj;
		for (i=0; i<histo->sizeX(); i++) for (j=0; j<histo->sizeY(); j++) histo->update(i,j,0.0f);
		for(i=0; i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++) if (mask.get(i,j))
		{
			bool cont=true;
				//printf("."); fflush(stdout);
			for (ii=0; ii<size && cont; ii++) for (jj=0; jj<size && cont; jj++)
			{
				int idx = i+ii-size/2, idy=j+jj-size/2;
				if (idx>=0 && idx<this->sizeX() && idy>=0 && idy<this->sizeY())
				{
					if (mask.get(idx,idy))
					{
						win.update(ii,jj,(float)this->get(idx,idy)/sc);
					}
					else cont=false;
				}
				else cont=false;
			}
			if (cont)
			{
				fft = win.fft(false, 1.0f, 0.0f);
				for (ii=0; ii<ff.size(); ii++)
				{
					hvVec2<int> freq=ff.at(ii);
					hvPair<float,float> vv = fft->get((freq.X()<0?size+freq.X():freq.X())+(freq.Y()<0?size+freq.Y():freq.Y())*size);
					float val = vv.getLeft()*vv.getLeft()+vv.getRight()*vv.getRight();
					if (freq.X()!=0 || freq.Y()!=0) val = (float)log(1.0+10000.0*val/0.1)/(float)log(10000.0+1.0);
					int ind = (int)(val*(float)histo->sizeX());
					if (ind<0) { std::cout<<"warning : out of bin range"<<val<<"\n"; ind=0; }
					else if (ind>=histo->sizeX()) { std::cout << "warning : out of bin range" << val << "\n"; ind=histo->sizeX()-1; }
					histo->update(ind,ii, histo->get(ind,ii)+1.0f);
				}
				nn++;
				delete fft;
			}
		}
		for (i=0; i<histo->sizeX(); i++) for (j=0; j<histo->sizeY(); j++) histo->update(i,j,5.0*histo->get(i,j)/(float)nn);
	}
	void wftpyramid(int nlevels, const hvBitmap &mask, int size, float sc, const std::vector<hvVec2<int> > &ff, hvPict<float> *histo[])
	{
		std::cout<<"level 0\n";
		wft(mask,size,sc,ff,histo[0]);
		hvBitmap mm,mm2; mm.shrink(mask,false);
		hvPict<T> pp,pp2; pp.shrink(this); 
		int i; for (i=1; i<nlevels; i++)
		{
			std::cout<<"level "<<i<<"\n";
			pp.wft(mm,size,sc,ff,histo[i]);
			mm2.shrink(mm,false); pp2.shrink(&pp);
			mm = mm2; pp.clone(pp2);
		}
	}
};


}

#endif // !efined(AFX_PICT_H__098453F0_1C38_49E9_A6F4_AABF90AA55E8__INCLUDED_)
