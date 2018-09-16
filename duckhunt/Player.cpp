#include <cstdlib>
#include <iostream>

#include "Player.hpp"
#include "matrix.h"
#include "globals.h"
#include "../HMM/hmm.h"
#include "classification.h"

namespace ducks
{

Player::Player()
{
    this->current_tstep = 0;
    this->current_round = 0;
    this->species_hmms = unordered_map<ESpecies, Lambda>(ESpecies::COUNT_SPECIES);
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
    }

    /*cerr << "Round\t" << pState.getRound()
        << "\tBirds\t" << pState.getNumBirds()
        << "\tMove\t" << pState.getBird(1).getLastObservation()
        << "\tNew turns\t" << pState.getNumNewTurns() << endl << flush;

    cerr << "current time: " << current_tstep << endl;*/
    size_t no_birds = pState.getNumBirds();

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

    int iters;
    //We wait some time before we start training our HMMs, to gather enough observations.
    if (this->current_tstep == 80) {
        for (int i = 0; i < no_birds; i++) {
            iters = hmm::model_estimate(this->HMMs[i], pDue);
            //cerr << "iterations bird " << i << ": " << iters << endl;
        }
    }

    if (current_tstep > 80)
        for (int i = 1; i < 2; i++)//i < pState.getNumBirds(); i++)
            if (pState.getBird(i).isAlive()) {
                /*cerr << this->HMMs[0].A << endl;
                cerr << this->HMMs[0].B << endl;
                cerr << this->HMMs[0].pi << endl;*/

                int guess = hmm::next_obs_guess(this->HMMs[i]);
                cerr << "Shooting at " << guess << endl;
                return Action(i, (EMovement) guess);
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
        cerr << "Guessing time!" << endl << flush;

        //note that the 0th element in the vector is always 0 as the actual groups are 1 <= n <= COUNT_SPECIES
        vector<int> optimal_guesses = classification::group_models(this->HMMs, ESpecies::COUNT_SPECIES, true);

        cerr << "Guesses calculated: ";

        for (int i : optimal_guesses)
            cerr << i << " ";
        
        cerr << endl;

        for (int i = 1; i < optimal_guesses.size(); i++)
            lGuesses[optimal_guesses[i]] = (ESpecies)(i - 1);
    }

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

    vector<bool> species_present(ESpecies::COUNT_SPECIES);

    cerr << "-- Actual species --" << endl;
    for (int i = 0; i < pSpecies.size(); i++)
        if (pSpecies[i] != ESpecies::SPECIES_UNKNOWN) {
            cerr << "Bird " << i << " was " << pSpecies[i] << endl;
            species_present[pSpecies[i]] = true;
        }
    
    int tot = 0;
    for (auto b : species_present)
        if (b)
            tot++;
    
    cerr << tot << " species in total" << endl;
}


} /*namespace ducks*/
