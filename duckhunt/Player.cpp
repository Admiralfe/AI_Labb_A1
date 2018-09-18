#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "Player.hpp"
#include "matrix.h"
#include "globals.h"
#include "hmm.h"
#include "classification.h"

namespace ducks
{

Player::Player()
{
    this->current_tstep = -1;
    this->current_round = -1;
    this->species_hmms = unordered_map<ESpecies, Lambda, std::hash<int>>(ESpecies::COUNT_SPECIES);
    this->shots = 0;
    this->hits = 0;
}

Action Player::shoot(const GameState &pState, const Deadline &pDue)
{
    /*
     * Here you should write your clever algorithms to get the best action.
     * This skeleton never shoots.
     */


    if (current_round != pState.getRound()) {
        current_round = pState.getRound();
        current_tstep = 0;
        if (current_round != 0)
            cerr << "Previous round took " << ((double) (clock() - roundTimer) / CLOCKS_PER_SEC * 1000.0) << "ms" << endl;
        roundTimer = clock();
        cerr << "Round " << current_round << endl;
        //cerr << "--Species HMMs--" << endl;

        shot_once = set<int>();
        shot_twice = set<int>();

        /*for (pair<ESpecies, Lambda> p : species_hmms) {
            cerr << "Observations for " << p.first << ": " << p.second.no_obs << "/" << p.second.obs_seq.size() << endl;
        }*/
    }

    /*cerr << "Round\t" << pState.getRound()
        << "\tBirds\t" << pState.getNumBirds()
        << "\tMove\t" << pState.getBird(1).getLastObservation()
        << "\tNew turns\t" << pState.getNumNewTurns() << endl << flush;

    cerr << "current time: " << current_tstep << endl;*/
    size_t no_birds = pState.getNumBirds();

    //cerr << "Birds: " << pState.getNumBirds() << endl;

    //Initialize an HMM for each bird on the first time step.
    if (this->current_tstep == 0) {
        this->HMMs = vector<Lambda>(no_birds);
        for (int i = 0; i < no_birds; i++) {
            this->HMMs[i] = Lambda();
        }
    }

    int no_new_turns = pState.getNumNewTurns();
    this->current_tstep += no_new_turns;

    //Add the new observations for each bird to the observation sequence.
    for (int i = 0; i < no_birds; i++) {
        for (int t = current_tstep - 1; t >= current_tstep - no_new_turns; t--) {
            if (pState.getBird(i).getObservation(t) != EMovement::MOVE_DEAD) {
                this->HMMs[i].obs_seq[t] = pState.getBird(i).getObservation(t);
                //this->HMMs[i].obs_seq.push_back(pState.getBird(i).getObservation(t));
                this->HMMs[i].no_obs++;
            }
        }
    }

    /*
    cerr << "observation sequence:" << endl;
    for (int t = 0; t < current_tstep; t++) {
        cerr << HMMs[1].obs_seq[t] << ", ";
    }

    cerr << endl;

    cerr << "number of observations: " << this->HMMs[0].no_obs << endl;
    */

    //cerr << "*";
    //cerr << "once: ";
    //for (int i : shot_once) cerr << i << " ";
    //cerr << endl;
    //cerr << "twice: ";
    //for (int i : shot_twice) cerr << i << " ";
    //cerr << endl;

    //We wait some time before we start training our HMMs, to gather enough observations.
    /*if (current_tstep == 100 - pState.getNumBirds()) {
        for (int i = 0; i < no_birds; i++) {
            hmm::model_estimate(this->HMMs[i], pDue);
            //cerr << "iterations bird " << i << ": " << iters << endl;
        }

    } else */if (current_tstep > 100 - pState.getNumBirds() 
                && pState.getRound() > 1
                && species_hmms.find(ESpecies::SPECIES_BLACK_STORK) != species_hmms.end()) { //Only want to train on first two rounds
        //cerr << "Estimating model parameters..." << endl;
        
        vector<tuple<ESpecies, number, int, number>> most_probable(no_birds);
        for (int i = 0; i < no_birds; i++) {
            most_probable[i] = { ESpecies::SPECIES_UNKNOWN, -numeric_limits<number>::infinity(), i, numeric_limits<number>::infinity()};

            if (pState.getBird(i).isAlive()) {
                hmm::model_estimate(this->HMMs[i], pDue, false, 20);
                vector<number> c(HMMs[i].no_obs);
                
                //cerr << "Probability for different species" << endl;
                for (int spec = 0; spec < ESpecies::COUNT_SPECIES; spec++) {
                    if (species_hmms.find((ESpecies) spec) == species_hmms.end())
                        continue;
                    Lambda mixed_model;
                    mixed_model.A = species_hmms[(ESpecies) spec].A;
                    mixed_model.B = species_hmms[(ESpecies) spec].B;
                    mixed_model.pi = species_hmms[(ESpecies) spec].pi;
                    mixed_model.obs_seq = HMMs[i].obs_seq;
                    mixed_model.no_obs = HMMs[i].no_obs;
                    for (int j = 0; j < c.size(); j++)
                        c[j] = 0;

                    hmm::a_pass(mixed_model, c);
                    number log_sum = 0;
                    
                    for (int j = 0; j < c.size(); j++)
                        log_sum -= log(c[j]);

                    //cerr << spec << ": " << log_sum << endl;
                    if ((ESpecies) spec == ESpecies::SPECIES_BLACK_STORK)
                        most_probable[i] = {
                            get<0>(most_probable[i]),
                            get<1>(most_probable[i]),
                            get<2>(most_probable[i]),
                            log_sum
                        };

                    if (!isnan(log_sum) && log_sum > get<1>(most_probable[i]))
                        most_probable[i] = {
                            (ESpecies) spec,
                            log_sum,
                            i,
                            get<3>(most_probable[i])
                        };
                }
            }
        }

        //cerr << "Guesses:" << endl;
        //for (auto p : most_probable)
        //    cerr << get<0>(p) << " ";
        //cerr << endl;

        sort(most_probable.begin(), most_probable.end(),
            [](tuple<ESpecies, number, int, number> a, tuple<ESpecies, number, int, number> b) {
                return isnan(get<3>(b)) || get<3>(a) < get<3>(b);
        });

        //cerr << "most_probable (sorted):" << endl;
        //for (auto p : most_probable)
        //    cerr << get<0>(p) << "\t" << get<1>(p) << "\t" << get<2>(p) << "\t" << get<3>(p) << "\t" << endl;

        int bird = -1;
        EMovement movement = EMovement::MOVE_DEAD;

        number a, b;

        for (int i = 0; i < most_probable.size(); i++) {
            auto tup = most_probable[i];
            if (pState.getBird(get<2>(tup)).isAlive()
                    && get<0>(tup) != ESpecies::SPECIES_BLACK_STORK
                    && get<0>(tup) != ESpecies::SPECIES_UNKNOWN
                    && shot_twice.find(get<2>(tup)) == shot_twice.end()
                    /*&& most_probable[i].second > något threshold*/) {
                //for (int j = 0; j < get<2>(tup); j++)
                //    cerr << "  ";
                //cerr << "^" << endl;
                
                //cerr << "Shooting at " << get<2>(tup) << ", which we believe is a " << get<0>(tup) << endl;
                bird = get<2>(tup);

                /*Lambda mixed_model;
                mixed_model.A = species_hmms[(ESpecies) most_probable[i].first].A;
                mixed_model.B = species_hmms[(ESpecies) most_probable[i].first].B;
                mixed_model.pi = species_hmms[(ESpecies) most_probable[i].first].pi;
                mixed_model.obs_seq = HMMs[i].obs_seq;
                mixed_model.no_obs = HMMs[i].no_obs;*/
                movement = (EMovement) hmm::next_obs_guess(HMMs[get<2>(tup)], a, b);
                if (b < 0.5)
                    continue;
                else
                    break;
            }
        }

        //Pick a bird that is alive and is the bird with the maximum probability to shoot.
        /*number log_prob = 0;
        number prob_sum = 0;
        number max_log_prob = -std::numeric_limits<number>::infinity();
        number prob = 0;
        number max_prob = 0;

        int bird = -1;
        int guess = -1;

        cerr << "Picking most likely bird" << endl;
        for (int i = 4; i < pState.getNumBirds(); i++) {
            int candidate_guess;
            if (pState.getBird(i).isAlive() && i != 10) {
                cerr << this->HMMs[0].A << endl;
                cerr << this->HMMs[0].B << endl;
                cerr << this->HMMs[0].pi << endl;
                candidate_guess = hmm::next_obs_guess(this->HMMs[i], log_prob, prob);
                cerr << "probability: " << prob << endl ;
                cerr << "log_prob : " << log_prob << endl << endl;
                if (prob > max_prob) {
                    guess = candidate_guess;
                    max_log_prob = log_prob;
                    bird = i;
                    max_prob = prob;
                }
            }
        }

        cerr << endl << "maximum probability: " << max_prob << endl;
        cerr << "maximum log probability: " << max_log_prob << endl;*/
        for (int i = 0; i < no_birds; i++) {
            this->HMMs[i].reset();
        }

        if (bird != -1 && b > 0.5) {
            if (shot_once.find(bird) != shot_once.end()) {
                shot_once.erase(bird);
                shot_twice.insert(bird);
            } else
                shot_once.insert(bird);

            //cerr << "Shooting at " << bird << " in direction " << movement << endl;
            cerr << "!";
            shots++;
            return Action(bird, movement);
        } else {
            cerr << "x";
            return cDontShoot;
        }
    } else {
        //train species models
        cerr << "*";
        if (current_round != 0) {
            ESpecies current = (ESpecies)(current_tstep % ESpecies::COUNT_SPECIES);
            if (species_hmms.find(current) != species_hmms.end())
                hmm::model_estimate(species_hmms[current], pDue, false, 10 - (current_round / 2));
        }
    }
    
    // This line choose not to shoot
    return cDontShoot;

    //This line would predict that bird 0 will move right and shoot at it
    //return Action(0, MOVE_RIGHT);
}

std::vector<ESpecies> Player::guess(const GameState &pState, const Deadline &pDue)
{
    /*
     * Here you should write your clever algorithms to guess the species of each bird.
     * This skeleton makes no guesses, better safe than sorry!
     */

    cerr << endl;

    std::vector<ESpecies> lGuesses(pState.getNumBirds(), SPECIES_UNKNOWN);

    if (current_round == 0)
        for (int i = 0; i < lGuesses.size(); i++)
            lGuesses[i] = (ESpecies)(i % ESpecies::COUNT_SPECIES);
    else {
        /*cerr << "Guessing time!" << endl << flush;

        //note that the 0th element in the vector is always 0 as the actual groups are 1 <= n <= COUNT_SPECIES
        vector<int> optimal_guesses = classification::group_models(this->HMMs, ESpecies::COUNT_SPECIES, true);

        cerr << "Guesses calculated: ";

        for (int i : optimal_guesses)
            cerr << i << " ";
        
        cerr << endl;

        for (int i = 1; i < optimal_guesses.size(); i++)
            lGuesses[optimal_guesses[i]] = (ESpecies)(i - 1);*/

        /*vector<ESpecies> undiscovered_species;
        for (int i = 0; i < ESpecies::COUNT_SPECIES; i++)
            if (species_hmms.find((ESpecies) i) == species_hmms.end())
                undiscovered_species.push_back((ESpecies) i);
        
        vector<pair<int, number>> total_distances(lGuesses.size());

        vector<int> reordering(HMMs[0].B.getHeight());

        for (int i = 0; i < lGuesses.size(); i++) {
            number tot = 0;
            number minimal = -1;
            ESpecies closest = ESpecies::SPECIES_UNKNOWN;

            for (int j = 0; j < ESpecies::COUNT_SPECIES; j++)
                if (species_hmms.find((ESpecies) j) != species_hmms.end()) {
                    //calculate distance in B and also the difference in state numbering
                    number distance = HMMs[i].B.distance_squared(species_hmms[(ESpecies) j].B, reordering, false);
                    //calculate distance in A, using the calculated state renumbering
                    distance += HMMs[i].A.distance_squared(species_hmms[(ESpecies) j].A, reordering, true);
                    
                    if (minimal == -1 || distance <= minimal) {
                        minimal = distance;
                        closest = (ESpecies) j;
                    }
                    tot += distance;
                }
            
            lGuesses[i] = closest;
            total_distances[i] = {i, tot};
        }

        sort(total_distances.begin(), total_distances.end(),
            [](pair<int, int> a, pair<int, int> b) {
                return a.second < b.second;
        });

        for (int i = 0; i < undiscovered_species.size() && lGuesses.size() - i - 1 >= 0; i++)
            lGuesses[i] = undiscovered_species[i];*/
        
        vector<pair<ESpecies, number>> most_probable(lGuesses.size());
        for (int i = 0; i < lGuesses.size(); i++) {
            hmm::model_estimate(this->HMMs[i], pDue);
            vector<number> c(HMMs[i].no_obs);

            most_probable[i] = { ESpecies::SPECIES_UNKNOWN, -numeric_limits<number>::infinity() };
            
            //cerr << "Probability for different species" << endl;
            for (int spec = 0; spec < ESpecies::COUNT_SPECIES; spec++) {
                if (species_hmms.find((ESpecies) spec) == species_hmms.end())
                    continue;
                Lambda mixed_model;
                mixed_model.A = species_hmms[(ESpecies) spec].A;
                mixed_model.B = species_hmms[(ESpecies) spec].B;
                mixed_model.pi = species_hmms[(ESpecies) spec].pi;
                mixed_model.obs_seq = HMMs[i].obs_seq;
                mixed_model.no_obs = HMMs[i].no_obs;

                for (int j = 0; j < c.size(); j++)
                    c[j] = 0;

                hmm::a_pass(mixed_model, c);
                number log_sum = 0;
                
                for (int j = 0; j < c.size(); j++)
                    log_sum -= log(c[j]);

                //cerr << spec << ": " << log_sum << endl;

                if (!isnan(log_sum) && log_sum > most_probable[i].second)
                    most_probable[i] = {(ESpecies) spec, log_sum};
            }
            lGuesses[i] = most_probable[i].first;
        }

    }

    cerr << "Guess: " << endl;
    for (auto i : lGuesses)
        cerr << i << " ";
    cerr << endl;

    return lGuesses;
}

void Player::hit(const GameState &pState, int pBird, const Deadline &pDue)
{
    /*
     * If you hit the bird you are trying to shoot, you will be notified through this function.
     */

    //std::cerr << "HIT BIRD!!!" << std::endl;
    hits++;
}

void Player::reveal(const GameState &pState, const std::vector<ESpecies> &pSpecies, const Deadline &pDue)
{
    /*
     * If you made any guesses, you will find out the true species of those birds in this function.
     */
    int storks = 0;
    cerr << "Actual species:" << endl;
    for (int i = 0; i < pSpecies.size(); i++)
        if (pSpecies[i] != ESpecies::SPECIES_UNKNOWN) {
            //cerr << "Bird " << i << " was " << pSpecies[i] << endl;
            cerr << pSpecies[i] << " ";

            if (pSpecies[i] == ESpecies::SPECIES_BLACK_STORK)
                storks++;

            if (species_hmms.find(pSpecies[i]) == species_hmms.end()) {
                //if we found a bird for which we have no model, make it the model of the species
                vector<int> trimmed(HMMs[i].no_obs);
                for (int j = 0; j < trimmed.size() && HMMs[i].obs_seq[j] != -1; j++)
                    trimmed[j] = HMMs[i].obs_seq[j];
                
                HMMs[i].obs_seq = trimmed;
                HMMs[i].no_obs = trimmed.size();
                species_hmms.insert(pair<ESpecies, Lambda>(pSpecies[i], HMMs[i]));
            } else {
                //else add its observations to the existing model
                int actualLength = 0;
                while (actualLength < HMMs[i].no_obs && HMMs[i].obs_seq[actualLength] != -1)
                    actualLength++;
                
                species_hmms[pSpecies[i]].obs_seq.insert(
                    species_hmms[pSpecies[i]].obs_seq.end(),
                    HMMs[i].obs_seq.begin(),
                    HMMs[i].obs_seq.begin() + actualLength);
                species_hmms[pSpecies[i]].no_obs += actualLength;

                species_hmms[pSpecies[i]].reset();
            }
        }
    cerr << endl;
    cerr << storks << " black stork" << (storks == 1 ? "" : "s") << endl;
    cerr << "So far hit " << hits << " of " << shots << " shot" << (shots == 1 ? "" : "s") << endl;
}


} /*namespace ducks*/
