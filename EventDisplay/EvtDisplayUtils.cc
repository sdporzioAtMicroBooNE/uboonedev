//
// - Serves as the link between ROOT "events" (e.g. mouse-clicks) and the ART
//   event display service by providing a receiver slot for signals generated
//   by the ROOT events.  A ROOT dictionary needs to be generated for this.
//

#include "EventDisplayBase/NavState.h"
#include "EvtDisplayUtils.h"
#include <string>

namespace evdb3D
{
    EvtDisplayUtils::EvtDisplayUtils()
        :fTbRun(0),fTbEvt(0)
    {}
  
    void EvtDisplayUtils::PrevEvent()
    {
        evdb::NavState::Set(evdb::kPREV_EVENT);
    }
    
    void EvtDisplayUtils::NextEvent()
    {
        evdb::NavState::Set(evdb::kNEXT_EVENT);
    }
    
    void EvtDisplayUtils::GotoEvent()
    {
        int run = std::stoi(fTbRun->GetString());
        int event = std::stoi(fTbEvt->GetString());
        evdb::NavState::SetTarget(run, event);
        evdb::NavState::Set(evdb::kGOTO_EVENT);
    }
    
}
