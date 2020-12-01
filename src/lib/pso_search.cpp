/*
        PSOVina version 2.0

        This file is revised from monte_carlo.cpp in AutoDock Vina.

        Authors: Giotto H. K. TAI  <giottotai@yahoo.com.hk>

                 Shirley W. I. SIU <shirleysiu@umac.mo>


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

#include "pso_search.h"
#include "coords.h"
#include "quasi_newton.h"
#include "pso_mutate.h"
#include "mutate.h"

output_type pso_search::operator()(model& m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator,int num_birds,int struggle, double w, double c1, double c2) const {
	output_container tmp;
	this->operator()(m, tmp, p, ig, p_widened, ig_widened, corner1, corner2, increment_me, generator, num_birds, struggle, w, c1, c2); // call the version that produces the whole container
	VINA_CHECK(!tmp.empty());
	return tmp.front();
}

bool metropolis_accept(fl old_f, fl new_f, fl temperature, rng& generator) {
	if(new_f < old_f) return true;
	const fl acceptance_probability = std::exp((old_f - new_f) / temperature);
	return random_fl(0, 1, generator) < acceptance_probability;
}

// out is sorted
void pso_search::operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator,int num_birds,int struggle, double w,double c1,double c2) const {
    vec authentic_v(1000, 1000, 1000); // FIXME? this is here to avoid max_fl/max_fl
    conf_size s = m.get_size();
	change g(s);
	output_type tmp(s, 0);
	tmp.c.randomize(corner1, corner2, generator);  //first randomize
    output_type newtmp(s, 0);
    newtmp.c.randomize(corner1, corner2, generator);  //first randomize
    fl best_e = max_fl;
	quasi_newton quasi_newton_par;
    quasi_newton_par.max_steps = ssd_par.evals;
	output_type tmp_rough = tmp;
	pso particle(num_birds, w, c1, c2,corner1,corner2,generator,tmp.c);
	double* PersonalBest = new double[100];
	for(int cou = 0; cou < 100; cou++)
		PersonalBest[cou] = 0;
	double energy = 0;
	int count = 0;
    int i;
    
	VINA_U_FOR(step, num_steps) {
		if(increment_me)
			++(*increment_me);
		output_type candidate = tmp_rough;
		output_type candidate_1 = tmp;

		pso_mutate_conf(candidate, candidate_1, m, mutation_amplitude, generator, &particle, PersonalBest, p, ig, g, hunt_cap, quasi_newton_par, step, min_rmsd, num_saved_mins, temperature, out, corner1, corner2); //for each particle loop
        tmp_rough = candidate;
        if(step == 0 || metropolis_accept(tmp.e, candidate_1.e, temperature, generator))
        {
            tmp = candidate_1;
            m.set(tmp.c); // FIXME? useless?
            if(tmp.e < best_e || out.size() < num_saved_mins) {
                quasi_newton_par(m, p, ig, tmp, g, authentic_v);
                m.set(tmp.c); // FIXME? useless?
                tmp.coords = m.get_heavy_atom_movable_coords();
                add_to_output_container(out, tmp, min_rmsd, num_saved_mins); // 20 - max size
                if(tmp.e < best_e)
                    best_e = tmp.e;
            }
        }
        if(std::abs(pso::gbest_fit - energy) > 0.0001)
        {
            output_type newtmp = candidate_1;
            for (i=0; i<3*particle.number; i++)
            {
                output_type tmp_markov = newtmp;
                markov_mutate_conf(tmp_markov, m, mutation_amplitude, generator, p, ig, g, hunt_cap, quasi_newton_par); //for each particle loop
                if(step == 0 || metropolis_accept(newtmp.e, tmp_markov.e, temperature, generator))
                {
                    newtmp = tmp_markov;
                    m.set(newtmp.c); // FIXME? useless?
                    if(newtmp.e < best_e || out.size() < num_saved_mins) {
                        quasi_newton_par(m, p, ig, newtmp, g, authentic_v);
                        m.set(newtmp.c); // FIXME? useless?
                        newtmp.coords = m.get_heavy_atom_movable_coords();
                        add_to_output_container(out, newtmp, min_rmsd, num_saved_mins); // 20 - max size
                        if(newtmp.e < best_e){
                            tmp = newtmp;
                            best_e = newtmp.e;
                            pso::gbest_fit = best_e;
                        }
                            
                    }
                }
            }
        }
        
        if(std::abs(pso::gbest_fit - energy) < 0.0001)
        {
            count += 1;
            if(count > struggle)
            { 
                step = num_steps; //break the loop
                count = 0;
            }
        }
        else
        {
			energy = pso::gbest_fit;
			count =0;
		}
	}
	VINA_CHECK(!out.empty());
	VINA_CHECK(out.front().e <= out.back().e); // make sure the sorting worked in the correct order
}
