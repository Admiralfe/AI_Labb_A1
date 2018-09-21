#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "Player.hpp"
#include "matrix.h"
#include "globals.h"
#include "hmm.h"

#define OBS_SEQ_SIZE 100
#define SAFETY_FACTOR 0.5
#define ROUNDS_NO_SHOOT 3
#define TRAINING_SESSIONS_PER_TURN 9

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
    this->start = clock();
}

Action Player::shoot(const GameState &pState, const Deadline &pDue)
{
    if (current_round != pState.getRound()) {
        current_round = pState.getRound();
        current_tstep = 0;
        if (current_round != 0)
            cerr << "Round took " << ((double) (clock() - roundTimer) / CLOCKS_PER_SEC * 1000.0) << "ms" << endl << endl;
        roundTimer = clock();
        cerr << "Round " << current_round << endl;

        shot_once = set<int>();
        shot_twice = set<int>();
        observations = vector<pair<vector<int>, int>>(pState.getNumBirds());
        for (int i = 0; i < observations.size(); i++)
            observations[i] = { vector<int>(OBS_SEQ_SIZE), 0 };
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

    /*vector containing bird species guesses
     * 0: estimated species,
     * 1: estimated next movement,
     * 2: bird number,
     * 3: logarithmic probability of having the estimated next movement
     */
    vector<tuple<ESpecies, EMovement, int, number>> most_probable(no_birds);

    if (current_tstep > 100 - (pState.getNumBirds() / SAFETY_FACTOR)
            //Only want to train on first few rounds
            && pState.getRound() >= ROUNDS_NO_SHOOT
            //Uncomment if you don't want to shoot if you haven't got a model for the stork
            //&& species_hmms.find(ESpecies::SPECIES_BLACK_STORK) != species_hmms.end()
        ) {
        for (int i = 0; i < no_birds; i++) {
            most_probable[i] = {
                    ESpecies::SPECIES_UNKNOWN,
                    EMovement::MOVE_DEAD,
                    i,
                    -numeric_limits<number>::infinity()
                };

            if (pState.getBird(i).isAlive()) {
                vector<number> c(observations[i].second);
                number prob_max = -numeric_limits<number>::infinity();
                
                for (int spec = 0; spec < ESpecies::COUNT_SPECIES; spec++) {
                    if (species_hmms.find((ESpecies) spec) == species_hmms.end())
                        continue;
                    
                    for (int j = 0; j < c.size(); j++)
                        c[j] = 0;
                    
                    //calculate probability of the bird belonging to species spec
                    delete hmm::a_pass(species_hmms[(ESpecies) spec], c, observations[i]);
                    number log_sum = 0;
                    
                    for (int j = 0; j < c.size(); j++)
                        log_sum += log(c[j]);
                    log_sum = -log_sum;

                    if (!isnan(log_sum) && log_sum > prob_max) {
                        prob_max = log_sum;
                        most_probable[i] = {
                            (ESpecies) spec,
                            get<1>(most_probable[i]),
                            i,
                            get<3>(most_probable[i])
                        };
                    }
                }

                number p;
                //calculate the next observation (i.e. movement) and its probability
                int move = hmm::next_obs_guess(species_hmms[get<0>(most_probable[i])], observations[get<2>(most_probable[i])], p);

                most_probable[i] = {
                        get<0>(most_probable[i]),
                        (EMovement) move,
                        get<2>(most_probable[i]),
                        p
                };
            }
        }

        //sort the list of guesses in increasing order of likelihood that the next observation guess is correct
        sort(most_probable.begin(), most_probable.end(),
            [](tuple<ESpecies, EMovement, int, number> a, tuple<ESpecies, EMovement, int, number> b) {
                return (isnan(get<3>(a)) && !isnan(get<3>(b))) || get<3>(a) < get<3>(b);
        });

        int bird = -1;
        EMovement movement = EMovement::MOVE_DEAD;

        number prob;

        for (int i = most_probable.size() - 1; i >= 0; i--) {
            auto tup = most_probable[i];
            //do not shoot at blacklisted species
            if (find(blacklist.begin(), blacklist.end(), get<0>(tup)) == blacklist.end()
                    && get<0>(tup) != ESpecies::SPECIES_UNKNOWN
                    //uncomment this line to shoot at most two times at a certain bird
                    //&& shot_twice.find(get<2>(tup)) == shot_twice.end()
                    ) {
                bird = get<2>(tup);
                movement = get<1>(tup);
                prob = get<3>(tup);

                if (prob < SAFETY_FACTOR)
                    continue;
                else
                    break;
            }
        }

        //only shoot if we are at least 55 % certain of hitting the bird
        if (bird != -1 && prob > SAFETY_FACTOR) {
            if (shot_once.find(bird) != shot_once.end()) {
                shot_once.erase(bird);
                shot_twice.insert(bird);
            } else
                shot_once.insert(bird);

            cerr << "|";
            if (current_tstep == 99)
                cerr << endl;

            shots++;
            return Action(bird, movement);
        } else {
            cerr << "x";
            if (current_tstep == 99)
                cerr << endl;
            return cDontShoot;
        }
    }
    else {
        //train species models
        cerr << "*";
        if (current_tstep == 99)
            cerr << endl;
        if (current_round != 0) {
            //train the model of one species per turn
            ESpecies current = (ESpecies)(current_tstep % ESpecies::COUNT_SPECIES);
            if (species_hmms.find(current) != species_hmms.end())
                hmm::model_estimate(
                    species_hmms[current],
                    species_observations[current],
                    false,
                    TRAINING_SESSIONS_PER_TURN - (current_round / 2)
                );
        }
    }
    
    return cDontShoot;
}

std::vector<ESpecies> Player::guess(const GameState &pState, const Deadline &pDue)
{
    bird_species_map.clear();

    std::vector<ESpecies> lGuesses(pState.getNumBirds(), SPECIES_UNKNOWN);

    if (current_round == 0) {
        //since pigeons are most prevalent, always guess pigeon for each bird the first round
        for (int i = 0; i < lGuesses.size(); i++)
            lGuesses[i] = ESpecies::SPECIES_PIGEON;

        blacklist = vector<ESpecies>();
        blacklist.push_back(SPECIES_BLACK_STORK);
    } else {
        //in essence the same guessing method as in shoot, but no movement prediction
        vector<pair<ESpecies, number>> most_probable(lGuesses.size());
        for (int i = 0; i < lGuesses.size(); i++) {
            vector<number> c(observations[i].second);

            most_probable[i] = { ESpecies::SPECIES_UNKNOWN, -numeric_limits<number>::infinity() };

            for (int spec = 0; spec < ESpecies::COUNT_SPECIES; spec++) {
                if (species_hmms.find((ESpecies) spec) == species_hmms.end())
                    continue;
                for (int j = 0; j < c.size(); j++)
                    c[j] = 0;

                delete hmm::a_pass(species_hmms[(ESpecies) spec], c, observations[i]);
                number log_sum = 0;

                for (int j = 0; j < c.size(); j++)
                    log_sum += log(c[j]);
                log_sum = -log_sum;

                if (!isnan(log_sum) && log_sum > most_probable[i].second)
                    most_probable[i] = {(ESpecies) spec, log_sum};
            }
            lGuesses[i] = most_probable[i].first;

            //save the guess for the reveal stage
            bird_species_map.insert({i, most_probable[i].first});
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
    hits++;
    if (current_tstep == 99)
        cerr << "\b\b!" << endl;
    else
        cerr << "\b!";
}

void Player::reveal(const GameState &pState, const std::vector<ESpecies> &pSpecies, const Deadline &pDue)
{
    int storks = 0;
    cerr << "Actual species:" << endl;
    for (int i = 0; i < pSpecies.size(); i++) {
        if (pSpecies[i] != ESpecies::SPECIES_UNKNOWN) {
            cerr << pSpecies[i] << " ";

            if (pSpecies[i] == ESpecies::SPECIES_BLACK_STORK) {
                storks++;

                //If the bird was a stork but we thought it was something else, add the species to a blacklist
                //since it is probably similar to the stork
                //However, since we are still very unsure in the first few rounds we don't add anything to blacklist then
                if (bird_species_map[i] != SPECIES_BLACK_STORK && pState.getRound() >= 2)
                    if (find(blacklist.begin(), blacklist.end(), bird_species_map[i]) == blacklist.end())
                        blacklist.push_back(bird_species_map[i]);
            }

            //add a species model if we've found a bird of that species
            if (species_hmms.find(pSpecies[i]) == species_hmms.end()) {
                species_hmms.insert({pSpecies[i], Lambda()});
                vector<pair<vector<int>, int>> observation_holder = vector<pair<vector<int>, int>>();
                species_observations.insert({pSpecies[i], observation_holder});
            }

            //add the observation sequence of the bird to the list of observations of the species
            species_observations[pSpecies[i]].push_back(observations[i]);

            //reset the species models so we can train them from scratch again
            species_hmms[pSpecies[i]].reset();
        }
    }
    cerr << endl;
    cerr << "Number of blacklisted species: " << blacklist.size() << endl;
    cerr << storks << " black stork" << (storks == 1 ? "" : "s") << endl;
    cerr << "So far hit " << hits << " of " << shots << " shot" << (shots == 1 ? "" : "s") << endl;
    if (pState.getRound() == 9)
        cerr << "Total time " << ((double) (clock() - start) / CLOCKS_PER_SEC * 1000.0) << "ms" << endl;
}


} /*namespace ducks*/
