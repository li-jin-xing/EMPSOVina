/*
        PSOVina version 2.0 

        Authors: Giotto H. K. TAI  <giottotai@yahoo.com.hk>

                 Shirley W. I. Siu <shirleysiu@umac.mo>


        Computational Biology and Bioinformatics Lab

        University of Macau

        http://cbbio.cis.umac.mo

*/
#include "pso.h"
#include "random.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

    // the vector of degree of freedom
	qt pso::gbest_orientation;
	fl* pso::gbest_torsion;
	vec pso::gbest_position;

	double pso::gbest_fit;

	pso::pso(int num_birds,double w,double c1,double c2,const vec corner1,const vec corner2, rng& g,conf& c)
	{
		  sz torsionSize = c.ligands[0].torsions.size();
		  this->w = w;
		  this->c1 = c1;
		  this->c2 = c2;
		  this->number = num_birds;
		  this->g = g;
		  this->corner1[0] = corner1[0];			//minmum
		  this->corner1[1] = corner1[1];
		  this->corner1[2] = corner1[2];
		  this->corner2[0] = corner2[0];			//maximum
		  this->corner2[1] = corner2[1];
		  this->corner2[2] = corner2[2];
		  this->torsionSize = (int)torsionSize;
		  this->R1Max_ = 1;
		  this->R1Min_ = 0;
		  this->R2Max_ = 1;
		  this->R2Min_ = 0;
		  pso::gbest_torsion = new fl[torsionSize];
		  init(g,c);
	}
	
	void pso::init(rng &g,conf& c)
	{

		int i;
		for(i=0;i<this->number;i++)
		{
			bird single_bird;
			single_bird.pbest_fit = 1.7976931348623158e+308;
			single_bird.tmp_fit = 1.7976931348623158e+308;
			single_bird.learning_rate = 1;
			single_bird.order = 0;
			//set position part
			single_bird.velocity = random_in_box(this->corner1,this->corner2,g);
			single_bird.current_position = random_in_box(this->corner1,this->corner2,g);
			
			//set orientation part
			single_bird.vO = random_inside_sphere(g);
			qt tmp_o = c.ligands[0].rigid.orientation;
			quaternion_increment(tmp_o,  random_inside_sphere(g));
			single_bird.current_orientation = tmp_o;
			
			//init. the array for the number of torsion
			single_bird.current_torsion=new fl[this->torsionSize];
			single_bird.vT=new fl[this->torsionSize];
			single_bird.pbest_torsion=new fl[this->torsionSize];
			
			for(int x=0;x<this->torsionSize;x++)						//init. all the torsion that the ligand has
			{
				single_bird.vT[x] = random_fl(-pi, pi, g);
				single_bird.current_torsion[x] = random_fl(-pi, pi, g);
			}
			
			particle.push_back(single_bird);					
		}
		
		
		pso::gbest_fit = 1.7976931348623158e+308;

	}
	
	void pso::computeQPSOPositions(int i, fl a)
	{
		//**********************************************************************************************************
		vec ave_position = this->calculateAveragePosition();
		for (int k=0; k<3; k++)
		{
			fl fi1, fi2, p, v, b, z;
			fi1 = random_fl(0,1,this->g);
			fi2 = random_fl(0,1,this->g);
			p = (fi1*particle[i].pbest_pos[k] + fi2*pso::gbest_position[k])/(fi1+fi2);
			v = std::log(1.0/random_fl(0,1,this->g));
			b = a * std::fabs(ave_position[k]-particle[i].current_position[k]);

			z = random_fl(0,1,this->g);
			if (z < 0.5)
				particle[i].current_position[k] = p + b * v;
			else
				particle[i].current_position[k] = p - b * v;
			
		}
		//***********************************************************************************************************
		//give a random position, if outside the search box
		if(particle[i].current_position[0] < corner1[0])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[1] < corner1[1])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[2] < corner1[2])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[0] > corner2[0])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[1] > corner2[1])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[2] > corner2[2])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
	}

    void pso::computeQPSOOrientation(int i, fl a)
	{
		//**********************************************************************************************************
		qt ave_orientation = this->calculateAverageOrientation();
		fl z;
		double d[4];
		qt fi1 = qt(random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g));
		qt fi2 = qt(random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g));
		qt p, v, b, u = qt(random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g));

		p = quaternion_div(quaternion_mul(fi1, particle[i].pbest_orientation)+quaternion_mul(fi2, pso::gbest_orientation), (fi1+fi2));
		v = quaternion_log(quaternion_div(qt(1.0,1.0,1.0,1.0), u));
		b = a * quaternion_abs(ave_orientation - particle[i].current_orientation);

		z = random_fl(0,1,this->g);
		if (z < 0.5)
			d[0] = p.R_component_1() + b.R_component_1() * v.R_component_1();
		else
			d[0] = p.R_component_1() - b.R_component_1() * v.R_component_1();
		
		z = random_fl(0,1,this->g);
		if (z < 0.5)
			d[1] = p.R_component_2() + b.R_component_2() * v.R_component_2();
		else
			d[1] = p.R_component_2() - b.R_component_2() * v.R_component_2();
		
		z = random_fl(0,1,this->g);
		if (z < 0.5)
			d[2] = p.R_component_3() + b.R_component_3() * v.R_component_3();
		else
			d[2] = p.R_component_3() - b.R_component_3() * v.R_component_3();
		
		z = random_fl(0,1,this->g);
		if (z < 0.5)
			d[3] = p.R_component_4() + b.R_component_4() * v.R_component_4();
		else
			d[3] = p.R_component_4() - b.R_component_4() * v.R_component_4();

		particle[i].current_orientation = qt(d[0], d[1], d[2], d[3]);
		quaternion_normalize_approx(particle[i].current_orientation);
		//***********************************************************************************************************
	}

	void pso::computeQPSOTorsion(int i ,sz which, fl a)
	{
		//**************************************************************************************************
		fl ave_torsion = this->calculateAverageTorsion(which);
		
		fl fi1, fi2, p, v, b, z;
		fi1 = random_fl(0,1,this->g);
		fi2 = random_fl(0,1,this->g);
		p = (fi1*particle[i].pbest_torsion[which] + fi2*pso::gbest_torsion[which])/(fi1+fi2);
		v = std::log(1.0/random_fl(0,1,this->g));
		b = a * std::fabs(ave_torsion-particle[i].current_torsion[which]);

		z = random_fl(0,1,this->g);
		if (z < 0.5)
			particle[i].current_torsion[which] = p + b * v;
		else
			particle[i].current_torsion[which] = p - b * v;

		//****************************************************************************************************
		
		// if(particle[i].current_torsion[which] > pi)
		// 	particle[i].current_torsion[which] = random_fl(-pi, pi, this->g);
		// else if(particle[i].current_torsion[which] < -pi)
		// 	particle[i].current_torsion[which] = random_fl(-pi, pi, this->g);
	}
	
	void pso::sortParticle(int num[],int a,int b)
	{
		double fitness[b-a+1];
		int i, j, tmp;
		for (i = 0; i < b-a+1; i++) {
			fitness[i] = this->getPersonalBest(i+a);
			num[i] = i+a;
		}

		for (i = 0; i < b-a+1 - 1; i++) {
			for (j = 1; j < b-a+1 - i; j++) {
				if (fitness[j - 1]< fitness[j]) {
					tmp = num[j - 1];
					num[j - 1] = num[j];
					num[j] = tmp;

					tmp = fitness[j - 1];
					fitness[j - 1] = fitness[j];
					fitness[j] = tmp;
				}
			}
		}
		for (i = 0; i < b-a+1; i++)
			particle[num[i]].order = i;
		for (i = a; i <= b; i++)
			particle[i].learning_rate = pow((1 - particle[i].order / double( b-a+1)), 0.01 * std::log( (b-a+1) / 8.0)); //M=8
	}

	void pso::computeSLPositions(int i, fl a, int num[], int start, int end)
	{

		int better;
		vec ave_position = this->calculateSLAveragePosition(start,end);
		for (int k=0; k<3; k++)
		{
			fl fi1, fi2, fi3, p, v, b, z;
			fi1 = random_fl(0,1,this->g);
			fi2 = random_fl(0,1,this->g);
			fi3 = random_fl(0,1,this->g);
			
			fl rand_learning = random_fl(0,1,this->g);
			if (rand_learning < particle[i].learning_rate) 
			{
				if (particle[i].order < end-start+1-1)
				better = num[rand() % (end-start+1 - particle[i].order -1 ) + particle[i].order +1];
				if (particle[i].order == end-start+1-1)
				better = num[end-start+1 - 1];
				p = (fi1 * particle[i].pbest_pos[k] + fi2 * fabs(particle[better].pbest_pos[k] - particle[i].pbest_pos[k]) + fi3 * fabs(ave_position[k] - particle[i].pbest_pos[k]))/(fi1 + fi2 + fi3);
			}

			v = std::log(1.0/random_fl(0,1,this->g));
			b = a * std::fabs(ave_position[k]-particle[i].current_position[k]);
			z = random_fl(0,1,this->g);
			if (z < 0.5)
				particle[i].current_position[k] = p + b * v;
			else
				particle[i].current_position[k] = p - b * v;
		}
		
		//give a random position, if outside the search box
		if(particle[i].current_position[0] < corner1[0])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[1] < corner1[1])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[2] < corner1[2])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[0] > corner2[0])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[1] > corner2[1])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[2] > corner2[2])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);

	}

    void pso::computeSLOrientation(int i, fl a, int num[], int start, int end)
	{
		//**********************************************************************************************************
		qt ave_orientation = this->calculateSLAverageOrientation(start, end);
		fl z;
		double d[4];
		qt fi1 = qt(random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g));
		qt fi2 = qt(random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g));
		qt fi3 = qt(random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g));
		qt p, v, b, u = qt(random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g));
		int better;
		fl rand_learning = random_fl(0,1,this->g);
		if (rand_learning < particle[i].learning_rate) 
		{
			if (particle[i].order < end-start+1 -1)
			better = num[rand() % (end-start+1 - particle[i].order -1 ) + particle[i].order +1];
			if (particle[i].order == end-start+1-1)
			better = num[end-start+1 - 1];
			p = quaternion_div(quaternion_mul(fi1, particle[i].pbest_orientation)+quaternion_mul(fi2, quaternion_abs(particle[better].pbest_orientation - particle[i].pbest_orientation))
							+quaternion_mul(fi3, quaternion_abs(ave_orientation -  particle[i].pbest_orientation)), (fi1 + fi2 + fi3));
		}
		v = quaternion_log(quaternion_div(qt(1.0,1.0,1.0,1.0), u));
		b = a * quaternion_abs(ave_orientation - particle[i].current_orientation);

		z = random_fl(0,1,this->g);
		if (z < 0.5)
			d[0] = p.R_component_1() + b.R_component_1() * v.R_component_1();
		else
			d[0] = p.R_component_1() - b.R_component_1() * v.R_component_1();
		
		z = random_fl(0,1,this->g);
		if (z < 0.5)
			d[1] = p.R_component_2() + b.R_component_2() * v.R_component_2();
		else
			d[1] = p.R_component_2() - b.R_component_2() * v.R_component_2();
		
		z = random_fl(0,1,this->g);
		if (z < 0.5)
			d[2] = p.R_component_3() + b.R_component_3() * v.R_component_3();
		else
			d[2] = p.R_component_3() - b.R_component_3() * v.R_component_3();
		
		z = random_fl(0,1,this->g);
		if (z < 0.5)
			d[3] = p.R_component_4() + b.R_component_4() * v.R_component_4();
		else
			d[3] = p.R_component_4() - b.R_component_4() * v.R_component_4();

		particle[i].current_orientation = qt(d[0], d[1], d[2], d[3]);
		quaternion_normalize_approx(particle[i].current_orientation);

	}


	void pso::computeSLTorsion(int i,sz which, fl a, int num[], int start, int end)
	{
		//**************************************************************************************************
		fl ave_torsion = this->calculateSLAverageTorsion(which, start, end);
		fl fi1, fi2, fi3, p, v, b, z;
		fi1 = random_fl(0,1,this->g);
		fi2 = random_fl(0,1,this->g);
		fi3 = random_fl(0,1,this->g);

		int better;
		fl rand_learning = random_fl(0,1,this->g);
		if (rand_learning < particle[i].learning_rate) 
		{
			if (particle[i].order < end-start+1 -1)
			better = num[rand() % (end-start+1 - particle[i].order -1 ) + particle[i].order +1];
			if (particle[i].order == end-start+1 -1)
			better = num[end-start+1 - 1];
			p = (fi1 * particle[i].pbest_torsion[which] + fi2 * fabs(particle[better].pbest_torsion[which] - particle[i].pbest_torsion[which]) + fi3 * fabs(ave_torsion - particle[i].pbest_torsion[which]))/(fi1 + fi2 + fi3);

		}
		v = std::log(1.0/random_fl(0,1,this->g));
		b = a * std::fabs(ave_torsion-particle[i].current_torsion[which]);

		z = random_fl(0,1,this->g);
		if (z < 0.5)
			particle[i].current_torsion[which] = p + b * v;
		else
			particle[i].current_torsion[which] = p - b * v;
	}
	
	void pso::computeRDPSOPositions(int i, fl alpha_cm, fl beta_cm)
	{
		//**********************************************************************************************************
		vec ave_position = this->calculateAveragePosition();
		fl fi1, fi2, p, v, z, rdn, V, alpha=0.75, beta=1;
		for (int k=0; k<3; k++)
		{
			fi1 = random_fl(0,1,this->g);
			fi2 = random_fl(0,1,this->g);
			p = (fi1*particle[i].pbest_pos[k] + fi2*pso::gbest_position[k])/(fi1+fi2);
			v = std::log(1.0/random_fl(0,1,this->g));
			z = random_fl(0,1,this->g);
			if (z < 0.5)
				rdn = v;
			else
				rdn = -v;

			V=alpha*fabs(ave_position[k]-particle[i].current_position[k])*rdn+beta*(p-particle[i].current_position[k]);

			if (V>(corner2[k]-corner1[k]))
				V=corner2[k]-corner1[k];
			if (V<(corner1[k]-corner2[k]))
				V=corner1[k]-corner2[k];
			particle[i].current_position[k] = particle[i].current_position[k] + V;
		}
		//***********************************************************************************************************

		//give a random position, if outside the search box
		if(particle[i].current_position[0] < corner1[0])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[1] < corner1[1])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[2] < corner1[2])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		
		if(particle[i].current_position[0] > corner2[0])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[1] > corner2[1])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[2] > corner2[2])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
	}

	void pso::computeRDPSOOrientation(int i, fl alpha_cm, fl beta_cm)
	{
		qt ave_orientation = this->calculateAverageOrientation();
		fl z, alpha=0.75, beta=1;
		double d[4];
		qt fi1 = qt(random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g));
		qt fi2 = qt(random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g));
		qt p, v, rdn, V, u = qt(random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g));
		p = quaternion_div(quaternion_mul(fi1, particle[i].pbest_orientation)+quaternion_mul(fi2, pso::gbest_orientation), (fi1+fi2));
		v = quaternion_log(quaternion_div(qt(1.0,1.0,1.0,1.0), u));
		z = random_fl(0,1,this->g);
		if (z < 0.5)
			rdn = v;
		else
			rdn = -v;

		V = alpha * quaternion_mul(quaternion_abs(ave_orientation - particle[i].current_orientation),rdn)
			+beta*(p-particle[i].current_orientation);
		particle[i].current_orientation = particle[i].current_orientation + V;
		quaternion_normalize_approx(particle[i].current_orientation);
	}

	void pso::computeRDPSOTorsion(int i,sz which, fl alpha_cm, fl beta_cm)
	{
		fl ave_torsion = this->calculateAverageTorsion(which);
		fl alpha=0.75, beta=1;
		fl fi1, fi2, p, rdn, v, z, V;
		fi1 = random_fl(0,1,this->g);
		fi2 = random_fl(0,1,this->g);
		p = (fi1*particle[i].pbest_torsion[which] + fi2*pso::gbest_torsion[which])/(fi1+fi2);
		v = std::log(1.0/random_fl(0,1,this->g));
		z = random_fl(0,1,this->g);
		if (z < 0.5)
			rdn = v;
		else
			rdn = -v;
		V=alpha*std::fabs(ave_torsion-particle[i].current_torsion[which])*rdn+beta*(p-particle[i].current_torsion[which]);
		particle[i].current_torsion[which] = particle[i].current_torsion[which] + V;
	}
	
	void pso::updateVelocity(int i,double cm,double l)
	{
			particle[i].velocity[0] = particle[i].velocity[0]*l+c1*cm*(particle[i].pbest_pos[0]-particle[i].current_position[0])+c2*(1-cm)*(pso::gbest_position[0]-particle[i].current_position[0]);
			particle[i].velocity[1] = particle[i].velocity[1]*l+c1*cm*(particle[i].pbest_pos[1]-particle[i].current_position[1])+c2*(1-cm)*(pso::gbest_position[1]-particle[i].current_position[1]);
			particle[i].velocity[2] = particle[i].velocity[2]*l+c1*cm*(particle[i].pbest_pos[2]-particle[i].current_position[2])+c2*(1-cm)*(pso::gbest_position[2]-particle[i].current_position[2]);
				
	}
	
	void pso::updateVelocityO(int i,double cm,double l)
	{
	
		    qt p1 = particle[i].pbest_orientation-particle[i].current_orientation;
		    qt p2 = pso::gbest_orientation-particle[i].current_orientation;
			particle[i].vO = particle[i].vO*l+c1*cm*quaternion_to_angle(p1)+c2*(1-cm)*quaternion_to_angle(p2);
			
	}
	
	void pso::updateVelocityT(int i,sz which,double cm,double l)
	{
			particle[i].vT[which] = particle[i].vT[which]*l+c1*cm*(particle[i].pbest_torsion[which]-particle[i].current_torsion[which])+c2*(1-cm)*(pso::gbest_torsion[which]-particle[i].current_torsion[which]);
	}

	void pso::computePSOPositions(int i)
	{

		particle[i].current_position[0] = particle[i].current_position[0] + particle[i].velocity[0];
		particle[i].current_position[1] = particle[i].current_position[1] + particle[i].velocity[1];
		particle[i].current_position[2] = particle[i].current_position[2] + particle[i].velocity[2];
		
		//give a random position, if outside the search box
		if(particle[i].current_position[0] < corner1[0])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[1] < corner1[1])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[2] < corner1[2])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		
		if(particle[i].current_position[0] > corner2[0])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[1] > corner2[1])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[2] > corner2[2])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);

	}

    void pso::computePSOOrientation(int i)
	{
			vec tmp_v = particle[i].vO;
			quaternion_increment(particle[i].current_orientation, tmp_v);
	}
	void pso::computePSOTorsion(int i, sz which)
	{
			particle[i].current_torsion[which] = particle[i].current_torsion[which] + particle[i].vT[which];

			//if(particle[i].current_torsion[which] > pi)
				//particle[i].current_torsion[which] = random_fl(-pi, pi, this->g);
			//else if(particle[i].current_torsion[which] < -pi)
				//particle[i].current_torsion[which] = random_fl(-pi, pi, this->g);
	}

	void pso::updatePersonalBest(int i,double e)
	{
		particle[i].pbest_fit = e;
	}
	
	void pso::updateGlobalBest(int i)
	{
		pso::gbest_fit = particle[i].pbest_fit;
	}
	
	void pso::updateGlobalBest_1(fl en)
	{
		pso::gbest_fit = en;
	}

	
	double pso::getPersonalBest(int i)
	{
		return particle[i].pbest_fit;
	}

	void pso::updateBestPosition(int i,vec pos)
	{
		particle[i].pbest_pos = pos;
		
	}
	
	void pso::updateBestOrientation(int i, qt orientation)
	{
		particle[i].pbest_orientation = orientation;
	}
	
	void pso::updateBestTorsion(int i, fl torsion,sz which)
	{
		particle[i].pbest_torsion[which] = torsion;
	}
//******************************************************************
	vec pso::calculateAveragePosition()
	{
		vec tmp = vec(0, 0, 0);
		for(int i=0;i<this->number;i++)
		{
			tmp[0] += particle[i].pbest_pos[0];
			tmp[1] += particle[i].pbest_pos[1];
			tmp[2] += particle[i].pbest_pos[2];
		}
		tmp[0] = tmp[0] / float(this->number);
		tmp[1] = tmp[1] / float(this->number);
		tmp[2] = tmp[2] / float(this->number);
		return tmp;
	}

	Eigen::MatrixXf dot(qt qutar) {
		Eigen::MatrixXf m(4, 4);
		float tmp[4] = { qutar.R_component_1(), qutar.R_component_2(), qutar.R_component_3(), qutar.R_component_4() };
		for (int i = 0; i < 4; i++) 
			for (int j = 0; j < 4; j++)
				m(i, j) = tmp[i] * tmp[j];
		
		return m;
	}

	qt pso::calculateAverageOrientation()
	{
		Eigen::MatrixXf M(4, 4), eig_vec;
		Eigen::MatrixXf eig_val;
		Eigen::MatrixXf::Index valMax;
		M << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		float f = 0;
		int len;
		len = this->number;
		for (int i = 0; i < len; i++)
			M += (1.0/len)*dot(particle[i].pbest_orientation);
		Eigen::EigenSolver<Eigen::MatrixXf> es(M);
		eig_val = es.pseudoEigenvalueMatrix();
		eig_vec = es.pseudoEigenvectors();
		eig_val.rowwise().sum().maxCoeff(&valMax);

		return qt(eig_vec(0, valMax),eig_vec(1, valMax),eig_vec(2, valMax),eig_vec(3,valMax));
	}
	
	fl pso::calculateAverageTorsion(sz which)
	{
		fl tmp = 0;
		for(int i=0; i<this->number; i++)
			tmp += particle[i].pbest_torsion[which];
		tmp = tmp / float(this->number);
		return tmp;
	}

	vec pso::calculateSLAveragePosition(int start, int end)
	{
		vec tmp = vec(0, 0, 0);
		for(int i=start;i<=end;i++)
		{
			tmp[0] += particle[i].pbest_pos[0];
			tmp[1] += particle[i].pbest_pos[1];
			tmp[2] += particle[i].pbest_pos[2];
		}
		tmp[0] = tmp[0] / float(end-start+1);
		tmp[1] = tmp[1] / float(end-start+1);
		tmp[2] = tmp[2] / float(end-start+1);
		return tmp;
	}

	qt pso::calculateSLAverageOrientation(int start, int end)
	{
		Eigen::MatrixXf M(4, 4), eig_vec;
		Eigen::MatrixXf eig_val;
		Eigen::MatrixXf::Index valMax;
		M << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		float f = 0;

		for (int i = start; i <= end; i++)
			M += (1.0/(end-start+1))*dot(particle[i].pbest_orientation);

		Eigen::EigenSolver<Eigen::MatrixXf> es(M);
		eig_val = es.pseudoEigenvalueMatrix();
		eig_vec = es.pseudoEigenvectors();
		eig_val.rowwise().sum().maxCoeff(&valMax);

		return qt(eig_vec(0, valMax),eig_vec(1, valMax),eig_vec(2, valMax),eig_vec(3,valMax));
	}
	

	fl pso::calculateSLAverageTorsion(sz which, int start, int end)
	{
		fl tmp = 0;
		for(int i=start; i<=end; i++)
			tmp += particle[i].pbest_torsion[which];
		tmp = tmp / float(end-start+1);
		return tmp;
	}

//******************************************************************
	vec pso::getCurrentPosition(int i)
	{
		return particle[i].current_position;
	}
	
	qt pso::getCurrentOrientation(int i)
	{
		return particle[i].current_orientation;
	}
	
	fl pso::getCurrentTorsion(int i,sz which)
	{
		return particle[i].current_torsion[which];
	}

