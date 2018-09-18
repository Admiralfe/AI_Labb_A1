#include <cstdlib>
#include <iostream>
#include <algorithm>

#include "Player.hpp"
#include "../HMM/matrix.h"
#include "../HMM/globals.h"
#include "../HMM/hmm.h"
#include "../HMM/classification.h"

namespace ducks
{

Player::Player()
{
    this->current_tstep = 0;
    this->current_round = 0;
    this->species_hmms = unordered_map<ESpecies, Lambda, std::hash<int>>(ESpecies::COUNT_SPECIES);
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
        cerr << "Round start" << endl;
        cerr << "--Species HMMs--" << endl;

        for (pair<ESpecies, Lambda> p : species_hmms) {
            cerr << "Observations for " << p.first << ": " << p.second.no_obs << "/" << p.second.obs_seq.size() << endl;
        }
    }

    /*cerr << "Round\t" << pState.getRound()
        << "\tBirds\t" << pState.getNumBirds()
        << "\tMove\t" << pState.getBird(1).getLastObservation()
        << "\tNew turns\t" << pState.getNumNewTurns() << endl << flush;

    cerr << "current time: " << current_tstep << endl;*/
    size_t no_birds = pState.getNumBirds();

    cerr << "Birds: " << pState.getNumBirds() << endl;

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
            this->HMMs[i].obs_seq[t] = pState.getBird(i).getObservation(t);
            //this->HMMs[i].obs_seq.push_back(pState.getBird(i).getObservation(t));
            this->HMMs[i].no_obs++;
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

    cerr << "time: " << current_tstep << endl;

    //We wait some time before we start training our HMMs, to gather enough observations.
    /*if (current_tstep == 100 - pState.getNumBirds()) {
        for (int i = 0; i < no_birds; i++) {
            hmm::model_estimate(this->HMMs[i], pDue);
            //cerr << "iterations bird " << i << ": " << iters << endl;
        }

    } else */if (current_tstep > 100 - pState.getNumBirds() && pState.getRound() != 0) { //Only want to train on first round.
        cerr << "Estimating model parameters..." << endl;
        for (int i = 0; i < no_birds; i++) {
            if (pState.getBird(i).isAlive()) {
                hmm::model_estimate(this->HMMs[i], pDue);
            }
        }

        //Pick a bird that is alive and is the bird with the maximum probability to shoot.
        number log_prob = 0;
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
                /*cerr << this->HMMs[0].A << endl;
                cerr << this->HMMs[0].B << endl;
                cerr << this->HMMs[0].pi << endl;*/
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
        cerr << "maximum log probability: " << max_log_prob << endl;

        for (int i = 0; i < no_birds; i++) {
            this->HMMs[i].reset();
        }

        if (max_prob > 0.5) {
            return Action(bird, (EMovement) guess);
        } else {
            return cDontShoot;
        }
    } else {
        //train species models
        if (current_round != 0) {
            ESpecies current = (ESpecies)(current_tstep % ESpecies::COUNT_SPECIES);
            if (species_hmms.find(current) != species_hmms.end())
                hmm::model_estimate(species_hmms[current], pDue);
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

        vector<ESpecies> undiscovered_species;
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
            lGuesses[i] = undiscovered_species[i];
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


    std::cerr << "HIT BIRD!!!" << std::endl;
}

void Player::reveal(const GameState &pState, const std::vector<ESpecies> &pSpecies, const Deadline &pDue)
{
    /*
     * If you made any guesses, you will find out the true species of those birds in this function.
     */

    cerr << "-- Actual species --" << endl;
    for (int i = 0; i < pSpecies.size(); i++)
        if (pSpecies[i] != ESpecies::SPECIES_UNKNOWN) {
            //cerr << "Bird " << i << " was " << pSpecies[i] << endl;
            cerr << pSpecies[i] << " ";

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
            }
        }
    
    //cerr << endl;
}


} /*namespace ducks*/
