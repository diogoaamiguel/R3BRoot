/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#ifndef R3BRpcPmtMappedData_H
#define R3BRpcPmtMappedData_H

#include "TObject.h"
#include <stdint.h>

class R3BRpcPmtMappedData : public TObject
{

  public:
    // Default Constructor
    R3BRpcPmtMappedData();

    /** Standard Constructor
     *@param channelId   Channel unique identifier
     *@param fineTime    Fine Time
     *@param coarseTime  Coarse Time
     *@param edge        Leading or Trailing
     *@param side        Side
     **/
    R3BRpcPmtMappedData(UShort_t channelId, uint64_t fineTime, uint64_t coarseTime, Int_t edge, UShort_t Side);

    // Destructor
    virtual ~R3BRpcPmtMappedData() {}

    // Getters
    inline UShort_t GetChannelId() const { return fChannelId; }
    inline uint64_t GetCoarseTime() const { return fCoarseTime; }
    inline uint64_t GetFineTime() const { return fFineTime; }
    inline Int_t GetEdge() const { return fEdge; }
    inline UShort_t GetSide() const { return fSide; }

  protected:
    UShort_t fChannelId;     // Channel unique identifier
    uint64_t fFineTime;      // Fine time
    uint64_t fCoarseTime;    // Coarse time
    Int_t fEdge;             // Leading or Trailing
    Short_t fSide;           // Top Or Bottom

  public:
    ClassDef(R3BRpcPmtMappedData, 1)
};

#endif