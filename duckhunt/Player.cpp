#include <cstdlib>
#include <iostream>

#include "Player.hpp"
#include "matrix.h"
#include "hmm.h"
#include "globals.h"

namespace ducks
{

Player::Player()
{

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
