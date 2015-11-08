#include "Target.h"
#include <iostream>

using namespace std;

// Initialize the total contacts static varible. Used to track the total number of times targets were found across Target objects.
int Target::total_contacts = 0;

Target::Target()
{
   // cout << "Target created" << endl;
  is_found = false;
  total_contacts = 0;
}

void Target::setFound()
{
    total_contacts++;
    is_found = true;
}

void Target::setNotFound()
{
    total_contacts = 0;
    is_found = false;
}

bool Target::isFound()
{
    return is_found;
}

int Target::getTotalContacts()
{
   // cout << total_contacts << endl;
    return total_contacts;
}
