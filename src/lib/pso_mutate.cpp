/*
       PSOVina version 2.0  

        This file is revised from monte_carlo.cpp in AutoDock Vina.

        Authors: Giotto H. K. TAI  <giottotai@yahoo.com.hk>

                 Shirley W. I. Siu <shirleysiu@umac.mo>


        Computational Biology and Bioinformatics Lab

        University of Macau

        http://cbbio.cis.umac.mo

*/
/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#include "pso_mutate.h"
#include "coords.h"
#include <math.h>
#include <time.h>
#define PI 3.14159265

bool pso_metropolis_accept(fl old_f, fl new_f, fl temperature, rng& generator) {
	if(new_f < old_f) return true;
	const fl acceptance_probability = std::exp((old_f - new_f) / temperature);
	return random_fl(0, 1, generator) < acceptance_probability;
}

sz pso_count_mutable_entities(const conf& c) {
	sz counter = 0;
	VINA_FOR_IN(i, c.ligands)
		counter += 2 + c.ligands[i].torsions.size();
	VINA_FOR_IN(i, c.flex)
		counter += c.flex[i].torsions.size();
	return counter;
}

void pso_mutate_conf(output_type& candidate, output_type& candidate_1, model& m, fl amplitude, rng& generator, pso* particle, double* PersonalBest, const precalculate& p ,const igrid& ig,change& g,const vec& v,quasi_newton& quasi_newton_par, int step, fl min_rmsd, sz num_saved_mins, fl temperature, output_container& out, const vec& corner1, const vec& corner2) { // ONE OF: 2A for position, similar amp for orientation, randomize torsion
	int shrink = 1; //rough factor:1/10 = 0.1
	output_type tmp_2 = candidate;
	output_type newtmp = tmp_2;
	output_type tmp_markov = newtmp;
	sz mutable_entities_num = pso_count_mutable_entities(candidate.c);
	if(mutable_entities_num == 0) return;
	int which_int = random_int(0, int(mutable_entities_num - 1), generator);
	VINA_CHECK(which_int >= 0);
	sz which = sz(which_int);
	VINA_CHECK(which < mutable_entities_num);
	
	srand((unsigned)time(NULL));
	double r_cm, w_cm, q_cm, sl_cm, alpha_cm, beta_cm;

	r_cm = (double) rand() / (RAND_MAX + 1.0);
	w_cm = (double) rand() / (RAND_MAX + 1.0);
	q_cm = (double) rand() / (RAND_MAX + 1.0);
	sl_cm = (double) rand() / (RAND_MAX + 1.0);
	alpha_cm = (double) rand() / (RAND_MAX + 1.0);
	beta_cm = (double) rand() / (RAND_MAX + 1.0);
	
	int num[particle->number];
	int y,a_floor,b_ceil;
	VINA_FOR_IN(i, candidate.c.ligands)
	{
		model tmp_m = m;
		const vec authentic_v(1000, 1000, 1000);
		//Take part the position (either position or orientation or torsion)
		if(which == 0) 
		{
			for (y=0;y<particle->number;y++)
			{
				candidate.c.ligands[i].rigid.position = particle->getCurrentPosition(y);
				candidate.c.ligands[i].rigid.orientation = particle->getCurrentOrientation(y);
				for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
					candidate.c.ligands[i].torsions[z] = particle->getCurrentTorsion(y,z);
				
				//rough local search
				quasi_newton_par(tmp_m, p, ig, candidate, g, v, shrink);
				
				//set the personal best(energy value and position);
				if(candidate.e < particle->getPersonalBest(y) || step <= 18) //Cr = 18
				{
					if(candidate.e < particle->getPersonalBest(y))
					{
						particle->updatePersonalBest(y,candidate.e);
						//***********************************************************************************************
						particle->updateBestPosition(y,candidate.c.ligands[i].rigid.position);
						particle->updateBestOrientation(y,candidate.c.ligands[i].rigid.orientation);
						for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
							particle->updateBestTorsion(y, candidate.c.ligands[i].torsions[z],z);
						//************************************************************************************************
					}
					tmp_2 = candidate;
					
					//full local search
					quasi_newton_par(tmp_m, p, ig, tmp_2, g, v);
					
					if (tmp_2.e < PersonalBest[y])
					{

						newtmp = tmp_2;
						for (int x=0; x<0.5*particle->number; x++)
						{	
							tmp_markov = newtmp;
							markov_mutate_conf(tmp_markov, tmp_m, amplitude, generator, p, ig, g, v, quasi_newton_par); //for each particle loop
							
							if(step == 0 || pso_metropolis_accept(newtmp.e, tmp_markov.e, temperature, generator))
							{
								newtmp = tmp_markov;
								tmp_m.set(newtmp.c); // FIXME? useless?
								if(newtmp.e < PersonalBest[y] || out.size() < num_saved_mins) {
									quasi_newton_par(tmp_m, p, ig, newtmp, g, authentic_v);
									tmp_m.set(newtmp.c); // FIXME? useless?
									newtmp.coords = tmp_m.get_heavy_atom_movable_coords();
									//m.set(newtmp.c);
									add_to_output_container(out, newtmp, min_rmsd, num_saved_mins); // 20 - max size
									if(newtmp.e < PersonalBest[y]){
										PersonalBest[y] = newtmp.e;
										tmp_2=newtmp;
									}
								}
							}
						}

						PersonalBest[y] = tmp_2.e;
						particle->updateBestPosition(y,tmp_2.c.ligands[i].rigid.position);
						particle->updateBestOrientation(y,tmp_2.c.ligands[i].rigid.orientation);
						for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
							particle->updateBestTorsion(y, tmp_2.c.ligands[i].torsions[z],z);
					}
						
					//set the global best(energy value and position);
					if(tmp_2.e < pso::gbest_fit)
					{
						particle->updateGlobalBest_1(tmp_2.e);
						
						pso::gbest_position = tmp_2.c.ligands[i].rigid.position;
						pso::gbest_orientation = tmp_2.c.ligands[i].rigid.orientation;
						for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
							pso::gbest_torsion[z] = tmp_2.c.ligands[i].torsions[z];
					}
				}

				if(y<particle->number*(2/12.0))
				{
					// q_cm = (double) rand() / (RAND_MAX + 1.0);
					//update chaotic map
					q_cm = 1.07*(7.86*q_cm-23.31*q_cm*q_cm+28.75*q_cm*q_cm*q_cm-13.302875*q_cm*q_cm*q_cm*q_cm);
					//compute the new Position;
					particle->computeQPSOPositions(y,q_cm);
				}
				else if(y<particle->number*(6/12.0))
				{
					// alpha_cm = (double) rand() / (RAND_MAX + 1.0);
					// beta_cm = (double) rand() / (RAND_MAX + 1.0);
					//compute the new Position;
					alpha_cm = 1.07*(7.86*alpha_cm-23.31*alpha_cm*alpha_cm+28.75*alpha_cm*alpha_cm*alpha_cm-13.302875*alpha_cm*alpha_cm*alpha_cm*alpha_cm);
					beta_cm = 1.07*(7.86*beta_cm-23.31*beta_cm*beta_cm+28.75*beta_cm*beta_cm*beta_cm-13.302875*beta_cm*beta_cm*beta_cm*beta_cm);
					particle->computeRDPSOPositions(y,alpha_cm,beta_cm);

				}
				else if(y<particle->number*(8/12.0))
				{
					// sl_cm = (double) rand() / (RAND_MAX + 1.0);
					a_floor=int(ceil(6/12.0*particle->number));
					b_ceil=int(floor(8/12.0*particle->number));
					sl_cm = 1.07*(7.86*sl_cm-23.31*sl_cm*sl_cm+28.75*sl_cm*sl_cm*sl_cm-13.302875*sl_cm*sl_cm*sl_cm*sl_cm);
					particle->sortParticle(num,a_floor,b_ceil);
					particle->computeSLPositions(y,sl_cm,num,a_floor,b_ceil);
				}
				else
				{	
					// r_cm = (double) rand() / (RAND_MAX + 1.0);
					// w_cm = (double) rand() / (RAND_MAX + 1.0);
					//update chaotic map
					r_cm = 1.07*(7.86*r_cm-23.31*r_cm*r_cm+28.75*r_cm*r_cm*r_cm-13.302875*r_cm*r_cm*r_cm*r_cm);
					w_cm = 1.07*(7.86*w_cm-23.31*w_cm*w_cm+28.75*w_cm*w_cm*w_cm-13.302875*w_cm*w_cm*w_cm*w_cm);
					
					//update each particle in every dimension
					particle->updateVelocity(y,r_cm,w_cm);
					
					//compute the new Position;
					particle->computePSOPositions(y);
				}

			}

			for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
				candidate_1.c.ligands[i].torsions[z] = pso::gbest_torsion[z];
			candidate_1.c.ligands[i].rigid.orientation = pso::gbest_orientation;
			candidate_1.c.ligands[i].rigid.position = pso::gbest_position;
			candidate_1.e = pso::gbest_fit;
			return; 
		}
			
			
		--which;
		//Take part orientation (either position or orientation or torsion)
		if(which == 0) 
		{
			fl gr = m.gyration_radius(i); 
			if(gr > epsilon_fl)
			{ // FIXME? just doing nothing for 0-radius molecules. do some other mutation?
		
				for (y=0;y<particle->number;y++)
				{
					candidate.c.ligands[i].rigid.position = particle->getCurrentPosition(y);
					candidate.c.ligands[i].rigid.orientation = particle->getCurrentOrientation(y);
					for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
						candidate.c.ligands[i].torsions[z] = particle->getCurrentTorsion(y,z);

					//rough local search
					quasi_newton_par(tmp_m, p, ig, candidate, g, v, shrink);

					//set the personal best(energy value and position);
					if(candidate.e < particle->getPersonalBest(y) || step <= 18)
					{
						if(candidate.e < particle->getPersonalBest(y))
						{
							particle->updatePersonalBest(y,candidate.e);
							//***********************************************************************************************
							particle->updateBestPosition(y,candidate.c.ligands[i].rigid.position);
							particle->updateBestOrientation(y,candidate.c.ligands[i].rigid.orientation);
							for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
								particle->updateBestTorsion(y, candidate.c.ligands[i].torsions[z],z);
							//***********************************************************************************************
						}
						tmp_2 = candidate;

						//full local search
						quasi_newton_par(tmp_m, p, ig, tmp_2, g, v);
						
						if(tmp_2.e < PersonalBest[y])
						{

							newtmp = tmp_2;
							for (int x=0; x<0.5*particle->number; x++)
							{	
								tmp_markov = newtmp;
								markov_mutate_conf(tmp_markov, tmp_m, amplitude, generator, p, ig, g, v, quasi_newton_par); //for each particle loop
								
								if(step == 0 || pso_metropolis_accept(newtmp.e, tmp_markov.e, temperature, generator))
								{
									newtmp = tmp_markov;
									tmp_m.set(newtmp.c); // FIXME? useless?
									if(newtmp.e < PersonalBest[y] || out.size() < num_saved_mins) {
										quasi_newton_par(tmp_m, p, ig, newtmp, g, authentic_v);
										tmp_m.set(newtmp.c); // FIXME? useless?
										newtmp.coords = tmp_m.get_heavy_atom_movable_coords();
										//m.set(newtmp.c);
										add_to_output_container(out, newtmp, min_rmsd, num_saved_mins); // 20 - max size
										if(newtmp.e < PersonalBest[y]){
											PersonalBest[y] = newtmp.e;
											tmp_2=newtmp;
										}
									}
								}
							}

							PersonalBest[y] = tmp_2.e;
							particle->updateBestPosition(y,tmp_2.c.ligands[i].rigid.position);
							particle->updateBestOrientation(y,tmp_2.c.ligands[i].rigid.orientation);
							for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
								particle->updateBestTorsion(y, tmp_2.c.ligands[i].torsions[z],z);
						}

						//set the global best(energy value and position);
						if(tmp_2.e< pso::gbest_fit)
						{
							particle->updateGlobalBest_1(tmp_2.e);						
							pso::gbest_position = tmp_2.c.ligands[i].rigid.position;
							pso::gbest_orientation = tmp_2.c.ligands[i].rigid.orientation;
							for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
								pso::gbest_torsion[z] = tmp_2.c.ligands[i].torsions[z];
						}
					}

					if(y<particle->number*(2/12.0))
					{
						// q_cm = (double) rand() / (RAND_MAX + 1.0);
						//update chaotic map
						q_cm = 1.07*(7.86*q_cm-23.31*q_cm*q_cm+28.75*q_cm*q_cm*q_cm-13.302875*q_cm*q_cm*q_cm*q_cm);
						//compute the new Orientation;
						particle->computeQPSOOrientation(y,q_cm);
					}
					else if(y<particle->number*(6/12.0))
					{
						// alpha_cm = (double) rand() / (RAND_MAX + 1.0);
						// beta_cm = (double) rand() / (RAND_MAX + 1.0);
						//compute the new Orientation;
						alpha_cm = 1.07*(7.86*alpha_cm-23.31*alpha_cm*alpha_cm+28.75*alpha_cm*alpha_cm*alpha_cm-13.302875*alpha_cm*alpha_cm*alpha_cm*alpha_cm);
						beta_cm = 1.07*(7.86*beta_cm-23.31*beta_cm*beta_cm+28.75*beta_cm*beta_cm*beta_cm-13.302875*beta_cm*beta_cm*beta_cm*beta_cm);
						particle->computeRDPSOOrientation(y,alpha_cm,beta_cm);

					}
					else if(y<particle->number*(8/12.0))
					{
						// sl_cm = (double) rand() / (RAND_MAX + 1.0);
						a_floor=int(ceil(6/12.0*particle->number));
						b_ceil=int(floor(8/12.0*particle->number));
						sl_cm = 1.07*(7.86*sl_cm-23.31*sl_cm*sl_cm+28.75*sl_cm*sl_cm*sl_cm-13.302875*sl_cm*sl_cm*sl_cm*sl_cm);
						particle->sortParticle(num,a_floor,b_ceil);
						particle->computeSLOrientation(y,sl_cm,num,a_floor,b_ceil);
					}
					else
					{
						// r_cm = (double) rand() / (RAND_MAX + 1.0);
						// w_cm = (double) rand() / (RAND_MAX + 1.0);
						//update chaotic map
						r_cm = 1.07*(7.86*r_cm-23.31*r_cm*r_cm+28.75*r_cm*r_cm*r_cm-13.302875*r_cm*r_cm*r_cm*r_cm);
						w_cm = 1.07*(7.86*w_cm-23.31*w_cm*w_cm+28.75*w_cm*w_cm*w_cm-13.302875*w_cm*w_cm*w_cm*w_cm);
						
						//update each particle in every dimension
						particle->updateVelocityO(y,r_cm,w_cm);
						//compute the new Orientation;
						particle->computePSOOrientation(y);
					}
				}

				for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
					candidate_1.c.ligands[i].torsions[z] = pso::gbest_torsion[z];
				candidate_1.c.ligands[i].rigid.orientation = pso::gbest_orientation;
				candidate_1.c.ligands[i].rigid.position = pso::gbest_position;
				candidate_1.e = pso::gbest_fit;
				return; 
			}
		}
		
		/*Torsions*/
		--which;
		if(which < candidate.c.ligands[i].torsions.size()) 
		{
			for (y=0;y<particle->number;y++)
			{
				candidate.c.ligands[i].rigid.position = particle->getCurrentPosition(y);
				candidate.c.ligands[i].rigid.orientation = particle->getCurrentOrientation(y);
				for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
					candidate.c.ligands[i].torsions[z] = particle->getCurrentTorsion(y,z);
				
				//rough local search
				quasi_newton_par(tmp_m, p, ig, candidate, g, v, shrink);

				//set the personal best(energy value and position);
				if(candidate.e < particle->getPersonalBest(y) || step <= 18)
				{
					if(candidate.e < particle->getPersonalBest(y))
					{
						particle->updatePersonalBest(y,candidate.e);
						//***********************************************************************************************
						particle->updateBestPosition(y,candidate.c.ligands[i].rigid.position);
						particle->updateBestOrientation(y,candidate.c.ligands[i].rigid.orientation);
						for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
							particle->updateBestTorsion(y, candidate.c.ligands[i].torsions[z],z);
						//***********************************************************************************************
					}
					tmp_2 = candidate;

					//full local search
					quasi_newton_par(tmp_m, p, ig, tmp_2, g, v);
					
					if(tmp_2.e < PersonalBest[y])
					{

						newtmp = tmp_2;
						for (int x=0; x<0.5*particle->number; x++)
						{	
							tmp_markov = newtmp;
							markov_mutate_conf(tmp_markov, tmp_m, amplitude, generator, p, ig, g, v, quasi_newton_par); //for each particle loop
							
							if(step == 0 || pso_metropolis_accept(newtmp.e, tmp_markov.e, temperature, generator))
							{
								newtmp = tmp_markov;
								tmp_m.set(newtmp.c); // FIXME? useless?
								if(newtmp.e < PersonalBest[y] || out.size() < num_saved_mins) {
									quasi_newton_par(tmp_m, p, ig, newtmp, g, authentic_v);
									tmp_m.set(newtmp.c); // FIXME? useless?
									newtmp.coords = tmp_m.get_heavy_atom_movable_coords();
									//m.set(newtmp.c);
									add_to_output_container(out, newtmp, min_rmsd, num_saved_mins); // 20 - max size
									if(newtmp.e < PersonalBest[y]){
										PersonalBest[y] = newtmp.e;
										tmp_2=newtmp;
									}
								}
							}
						}

						PersonalBest[y] = tmp_2.e;
						particle->updateBestPosition(y,tmp_2.c.ligands[i].rigid.position);
						particle->updateBestOrientation(y,tmp_2.c.ligands[i].rigid.orientation);
						for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
							particle->updateBestTorsion(y, tmp_2.c.ligands[i].torsions[z],z);
					}
					
					
					//set the global best(energy value and position);
					if(tmp_2.e < pso::gbest_fit)
					{
						particle->updateGlobalBest_1(tmp_2.e);
						pso::gbest_position = tmp_2.c.ligands[i].rigid.position;
						pso::gbest_orientation = tmp_2.c.ligands[i].rigid.orientation;
						for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
							pso::gbest_torsion[z] = tmp_2.c.ligands[i].torsions[z];
					}
				}

				if(y<particle->number*(2/12.0))
				{
					// q_cm = (double) rand() / (RAND_MAX + 1.0);
					//update chaotic map
					q_cm = 1.07*(7.86*q_cm-23.31*q_cm*q_cm+28.75*q_cm*q_cm*q_cm-13.302875*q_cm*q_cm*q_cm*q_cm);
					//compute the new Torsion;
					particle->computeQPSOTorsion(y,which,q_cm);
				}
				else if(y<particle->number*(6/12.0))
				{
					// alpha_cm = (double) rand() / (RAND_MAX + 1.0);
					// beta_cm = (double) rand() / (RAND_MAX + 1.0);
					//compute the new Torsion;
					alpha_cm = 1.07*(7.86*alpha_cm-23.31*alpha_cm*alpha_cm+28.75*alpha_cm*alpha_cm*alpha_cm-13.302875*alpha_cm*alpha_cm*alpha_cm*alpha_cm);
					beta_cm = 1.07*(7.86*beta_cm-23.31*beta_cm*beta_cm+28.75*beta_cm*beta_cm*beta_cm-13.302875*beta_cm*beta_cm*beta_cm*beta_cm);
					particle->computeRDPSOTorsion(y,which,alpha_cm,beta_cm);
				}
				else if(y<particle->number*(8/12.0))
				{
					// sl_cm = (double) rand() / (RAND_MAX + 1.0);
					a_floor=int(ceil(6/12.0*particle->number));
					b_ceil=int(floor(8/12.0*particle->number));
					particle->sortParticle(num,a_floor,b_ceil);
					sl_cm = 1.07*(7.86*sl_cm-23.31*sl_cm*sl_cm+28.75*sl_cm*sl_cm*sl_cm-13.302875*sl_cm*sl_cm*sl_cm*sl_cm);
					particle->computeSLTorsion(y,which,sl_cm,num,a_floor,b_ceil);
				}
				else
				{
					// r_cm = (double) rand() / (RAND_MAX + 1.0);
					// w_cm = (double) rand() / (RAND_MAX + 1.0);
					//update chaotic map
					r_cm = 1.07*(7.86*r_cm-23.31*r_cm*r_cm+28.75*r_cm*r_cm*r_cm-13.302875*r_cm*r_cm*r_cm*r_cm);
					w_cm = 1.07*(7.86*w_cm-23.31*w_cm*w_cm+28.75*w_cm*w_cm*w_cm-13.302875*w_cm*w_cm*w_cm*w_cm);
					//update each particle in every dimension
					particle->updateVelocityT(y,which,r_cm,w_cm);
					//compute the new Torsion;
					particle->computePSOTorsion(y,which);
				}

			}
			
			for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
				candidate_1.c.ligands[i].torsions[z] = pso::gbest_torsion[z];
			candidate_1.c.ligands[i].rigid.orientation = pso::gbest_orientation;
			candidate_1.c.ligands[i].rigid.position = pso::gbest_position;
			candidate_1.e = pso::gbest_fit;
			return; 
		}
		which -= candidate.c.ligands[i].torsions.size();
	}
		
	VINA_FOR_IN(i, candidate.c.flex) {
		if(which < candidate.c.flex[i].torsions.size()) { candidate.c.flex[i].torsions[which] = random_fl(-pi, pi, generator); return; }
		which -= candidate.c.flex[i].torsions.size();
	}
}

void markov_mutate_conf(output_type& candidate, const model& m, fl amplitude, rng& generator, const precalculate& p ,const igrid& ig,change& g,const vec& v,quasi_newton& quasi_newton_par)
{ // ONE OF: 2A for position, similar amp for orientation, randomize torsion
	
	sz mutable_entities_num = pso_count_mutable_entities(candidate.c);
	if(mutable_entities_num == 0) return;
	int which_int = random_int(0, int(mutable_entities_num - 1), generator);
	VINA_CHECK(which_int >= 0);
	sz which = sz(which_int);
	VINA_CHECK(which < mutable_entities_num);

	VINA_FOR_IN(i, candidate.c.ligands)
	{
		model tmp_m = m;
		if(which == 0) 
			candidate.c.ligands[i].rigid.position += amplitude * random_inside_sphere(generator); 
		--which;

		if(which == 0)
		{
			fl gr = m.gyration_radius(i); 
			if(gr > epsilon_fl) { // FIXME? just doing nothing for 0-radius molecules. do some other mutation?
				vec rotation; 
				rotation = amplitude / gr * random_inside_sphere(generator); 
				quaternion_increment(candidate.c.ligands[i].rigid.orientation, rotation);
			}
		}
		--which;

		if(which < candidate.c.ligands[i].torsions.size()) 
		{
			candidate.c.ligands[i].torsions[which] = random_fl(-pi, pi, generator); 
		}
		which -= candidate.c.ligands[i].torsions.size();

		quasi_newton_par(tmp_m, p, ig, candidate, g, v);

	}

	VINA_FOR_IN(i, candidate.c.flex) {
		if(which < candidate.c.flex[i].torsions.size()) { candidate.c.flex[i].torsions[which] = random_fl(-pi, pi, generator); return; }
		which -= candidate.c.flex[i].torsions.size();
	}
}