//#define FELIX_EFTERBLIVNA_CLION

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <cmath>

#ifdef FELIX_EFTERBLIVNA_CLION
#include "Player.hpp"
#include "../HMM/matrix.h"
#include "../HMM/globals.h"
#include "../HMM/hmm.h"
#include "../HMM/classification.h"
#else
#include "Player.hpp"
#include "matrix.h"
#include "globals.h"
#include "hmm.h"
#include "classification.h"
#endif

#define OBS_SEQ_INITIAL_SIZE 1000

namespace ducks
{

Player::Player()
{
    this->current_tstep = -1;
    this->current_round = -1;
    this->species_hmms = unordered_map<ESpecies, Lambda, std::hash<int>>(ESpecies::COUNT_SPECIES);
    this->species_total_observations = unordered_map<ESpecies, pair<vector<int>, int>, std::hash<int>>(ESpecies::COUNT_SPECIES);
    this->shots = 0;
    this->hits = 0;
}

Action Player::shoot(const GameState &pState, const Deadline &pDue)
{
    if (current_round != pState.getRound()) {
        current_round = pState.getRound();
        current_tstep = 0;
        if (current_round != 0)
            cerr << "Previous round took " << ((double) (clock() - roundTimer) / CLOCKS_PER_SEC * 1000.0) << "ms" << endl;
        roundTimer = clock();
        cerr << "Round " << current_round << endl;

        shot_once = set<int>();
        shot_twice = set<int>();
        observations = vector<pair<vector<int>, int>>(pState.getNumBirds());
        for (int i = 0; i < observations.size(); i++)
            observations[i] = { vector<int>(100), 0 };
    }

    size_t no_birds = pState.getNumBirds();

    int no_new_turns = pState.getNumNewTurns();
    this->current_tstep += no_new_turns;

    //Add the new observations for each bird to the observation sequence.
    for (int i = 0; i < no_birds; i++) {
        for (int t = current_tstep - 1; t >= current_tstep - no_new_turns; t--) {
            if (pState.getBird(i).getObservation(t) != EMovement::MOVE_DEAD) {
                observations[i].first[t] = pState.getBird(i).getObservation(t);
                observations[i].second++;
            }
        }
    }
    
    if (current_tstep > 100 - (pState.getNumBirds() * 1.5)
                //Only want to train on first two rounds
                && pState.getRound() > 1
                && species_hmms.find(ESpecies::SPECIES_BLACK_STORK) != species_hmms.end()) {
        //cerr << "Estimating model parameters..." << endl;
        
        vector<tuple<ESpecies, number, int, number>> most_probable(no_birds);
        for (int i = 0; i < no_birds; i++) {
            most_probable[i] = { ESpecies::SPECIES_UNKNOWN, -numeric_limits<number>::infinity(), i, numeric_limits<number>::infinity()};

            if (pState.getBird(i).isAlive()) {
                vector<number> c(observations[i].second);
                
                //cerr << "Probability for different species" << endl;
                for (int spec = 0; spec < ESpecies::COUNT_SPECIES; spec++) {
                    if (species_hmms.find((ESpecies) spec) == species_hmms.end())
                        continue;
                    for (int j = 0; j < c.size(); j++)
                        c[j] = 0;

                    hmm::a_pass(species_hmms[(ESpecies) spec], c, observations[i]);
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

        sort(most_probable.begin(), most_probable.end(),
            [](tuple<ESpecies, number, int, number> a, tuple<ESpecies, number, int, number> b) {
                return isnan(get<3>(b)) || get<3>(a) < get<3>(b);
        });

        int bird = -1;
        EMovement movement = EMovement::MOVE_DEAD;

        number a, prob;

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
                movement = (EMovement) hmm::next_obs_guess(species_hmms[get<0>(tup)], observations[get<2>(tup)], prob);
                if (prob < 0.5)
                    continue;
                else
                    break;
            }
        }

        if (bird != -1 && prob > 0.5) {
            if (shot_once.find(bird) != shot_once.end()) {
                shot_once.erase(bird);
                shot_twice.insert(bird);
            } else
                shot_once.insert(bird);

            //cerr << "Shooting at " << bird << " in direction " << movement << endl;
            cerr << "|";
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
                hmm::model_estimate(species_hmms[current], false, 10 - (current_round / 2));
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
            vector<number> c(observations[i].second);

            most_probable[i] = { ESpecies::SPECIES_UNKNOWN, -numeric_limits<number>::infinity() };
            
            for (int spec = 0; spec < ESpecies::COUNT_SPECIES; spec++) {
                if (species_hmms.find((ESpecies) spec) == species_hmms.end())
                    continue;
                for (int j = 0; j < c.size(); j++)
                    c[j] = 0;

                hmm::a_pass(species_hmms[(ESpecies) spec], c, observations[i]);
                number log_sum = 0;
                
                for (int j = 0; j < c.size(); j++)
                    log_sum -= log(c[j]);

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
    cerr << "\b!";
}

void Player::reveal(const GameState &pState, const std::vector<ESpecies> &pSpecies, const Deadline &pDue)
{
    int storks = 0;
    cerr << "Actual species:" << endl;
    for (int i = 0; i < pSpecies.size(); i++)
        if (pSpecies[i] != ESpecies::SPECIES_UNKNOWN) {
            //cerr << "Bird " << i << " was " << pSpecies[i] << endl;
            cerr << pSpecies[i] << " ";

            if (pSpecies[i] == ESpecies::SPECIES_BLACK_STORK)
                storks++;

            if (species_hmms.find(pSpecies[i]) == species_hmms.end()) {
                species_hmms.insert({pSpecies[i], Lambda()});
                species_total_observations.insert({pSpecies[i], {vector<int>(OBS_SEQ_INITIAL_SIZE), 0}});
            
                //else add its observations to the existing model
                int actualLength = 0;
                while (actualLength < observations[i].second && observations[i].first[actualLength] != -1)
                    actualLength++;
                
                species_total_observations[pSpecies[i]].first.insert(
                    species_total_observations[pSpecies[i]].first.end(),
                    observations[i].first.begin(),
                    observations[i].first.begin() + actualLength);
                species_total_observations[pSpecies[i]].second += actualLength;

                species_hmms[pSpecies[i]].reset();
                
            }
        }
    cerr << endl;
    cerr << storks << " black stork" << (storks == 1 ? "" : "s") << endl;
    cerr << "So far hit " << hits << " of " << shots << " shot" << (shots == 1 ? "" : "s") << endl;
}


} /*namespace ducks*/
