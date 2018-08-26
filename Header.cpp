#include "Header.h"

/* Internal parameters
- PARAM_INTERNAL_SAVE: save to disk
*/
bool PARAM_INTERNAL_SAVE;


/**
Driver function to ascent sort the vector of pairs by the second element of pairs
**/
bool ascentSortBySec(const IDPair &a, const IDPair &b)
{
    return (a.second < b.second);
}

/**
Driver function to descent sort the vector of pairs by the second element of pairs
**/
bool descentSortBySec(const IDPair &a, const IDPair &b)
{
    return (a.second > b.second);
}
