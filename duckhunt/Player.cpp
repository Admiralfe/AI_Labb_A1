#include <cstdlib>
#include <iostream>

#include "Player.hpp"
#include "matrix.h"
#include "globals.h"
#include "../HMM/hmm.h"

namespace ducks
{

Player::Player()
{
    this->current_tstep = 0;
}

Action Player::shoot(const GameState &pState, const Deadline &pDue)
{
    /*
     * Here you should write your clever algorithms to get the best action.
     * This skeleton never shoots.
     */

    cerr << "Round\t" << pState.getRound()
        << "\tBirds\t" << pState.getNumBirds()
        << "\tMove\t" << pState.getBird(0).getLastObservation()
        << "\tNew turns\t" << pState.getNumNewTurns() << endl << flush;

    int no_new_turns = pState.getNumNewTurns();
    this->current_tstep += no_new_turns;
    size_t no_birds = pState.getNumBirds();

    //Initialize an HMM for each bird on the first time step.
    if (this->current_tstep == 0) {
        for (int i = 0; i < no_birds; i++) {
            this->HMMs[i] = init_lambda();
        }
    }

    //Add the new observations for each bird to the observation sequence.
    for (int i = 0; i < no_birds; i++) {
        for (int t = current_tstep; t >= current_tstep - no_new_turns; t--) {
            this->HMMs[i].obs_seq[t] = pState.getBird(i).getObservation(t);
        }
    }

    //We wait some time before we start training our HMMs, to gather enough observations.
    if (this->current_tstep >= 14) {
        for (int i = 0; i < no_birds; i++) {
            hmm::model_estimate(this->HMMs[i];
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

    cerr << "Guessing time!" << endl << flush;

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
}


} /*namespace ducks*/
