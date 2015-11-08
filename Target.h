#ifndef TARGET_H
#define TARGET_H

#include "Agent.h"

class Target : public Agent
{
public:
    Target();
    bool isFound();
    void setFound();
    void setNotFound();
    int getTotalContacts();


private:
    static int total_contacts;
    bool is_found;
};

#endif // TARGET_H
