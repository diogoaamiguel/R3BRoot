/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2023 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "R3BDigitizingTamex.h"

#include <FairRuntimeDb.h>
#include <cmath>

#include "R3BNeulandHitModulePar.h"
#include "R3BNeulandHitPar.h"
#include <FairRunAna.h>
#include <R3BLogger.h>
#include <algorithm>
#include <iostream>
#include <utility>

namespace R3B::Digitizing::Neuland::Tamex
{
    // some declarations for static functions:
    namespace
    {
        template <class T>
        auto CheckOverlapping(const T& peak, std::vector<T>& peaks) -> decltype(peaks.begin());
        template <class T>
        void ReOverlapping(typename std::vector<T>::iterator v_iter, std::vector<T>& peaks);
        template <class T>
        void RemovePeakAt(typename std::vector<T>::iterator v_iter, std::vector<T>& peaks);
    } // namespace

    // global variables for default options:
    const size_t TmxPeaksInitialCapacity = 10;
    const double PMTPeak::peakWidth = 15.0;              // ns
    R3BNeulandHitPar* Channel::fNeulandHitPar = nullptr; // NOLINT

    Params::Params(TRandom3& rnd)
        : fRnd{ &rnd }
    {
    }

    Params::Params(const Params& other)
    {
        if (other.fRnd == nullptr)
        {
            throw std::runtime_error(
                "R3BDigitizingTamex: copy constructor of Params cannot take nullptr of random generator!");
        }
        *this = other;
    }

    PMTPeak::PMTPeak(Digitizing::Channel::Hit pmtHit, const Channel& channel)
    {
        auto par = channel.GetParConstRef();
        // apply saturation coefficent
        if (par.fExperimentalDataIsCorrectedForSaturation)
        {
            fQdc = pmtHit.light / (1. + par.fSaturationCoefficient * pmtHit.light);
        }
        else
        {
            fQdc = pmtHit.light;
        }
        fLETime = pmtHit.time;
    }

    auto PMTPeak::operator+=(const PMTPeak& rhs) -> PMTPeak&
    {
        fQdc += rhs.fQdc;
        fLETime = (fLETime < rhs.fLETime) ? fLETime : rhs.fLETime;
        return *this;
    }

    Peak::Peak(const PMTPeak& pmtPeak, Channel* channel)
        : fQdc(pmtPeak.GetQDC())
        , fLETime(pmtPeak.GetLETime())
        , fChannel(channel)
    {
        if (fChannel == nullptr)
        {
            LOG(fatal) << "channel is not bound to FQTPeak object!";
        }
        const auto& par = channel->GetParConstRef();

        // calculate the time and the width of the signal
        fWidth = QdcToWidth(fQdc, par);
        fTETime = fLETime + fWidth;
    }

    auto Peak::operator==(const Peak& sig) const -> bool
    {
        if (sig.fLETime == 0 && fLETime == 0)
        {
            LOG(warn) << "the times of both PMT signals are 0!";
        }
        return (fLETime <= (sig.fLETime + sig.fWidth)) && (sig.fLETime <= (fLETime + fWidth));
    }

    void Peak::operator+=(const Peak& sig)
    {
        if (fChannel == nullptr)
        {
            LOG(fatal) << "channel is not bound to FQTPeak object!";
        }
        fLETime = (fLETime < sig.fLETime) ? fLETime : sig.fLETime;
        fTETime = (fTETime > sig.fTETime) ? fTETime : sig.fTETime;
        fWidth = fTETime - fLETime;
        fQdc = WidthToQdc(fWidth, fChannel->GetParConstRef());
    }

    Channel::Channel(ChannelSide side, TRandom3& rnd)
        : Channel{ side, Params{ rnd } }
    {
    }

    Channel::Channel(ChannelSide side, const Params& par)
        : Digitizing::Channel{ side }
        , par_{ par }
    {
        fPMTPeaks.reserve(TmxPeaksInitialCapacity);
    }

    void Channel::GetHitPar(const std::string& hitParName)
    {
        if (hitParName.empty())
        {
            LOG(info) << "DigitizingTamex: Using default parameter for Tamex Channels.";
            return;
        }
        auto* run = FairRunAna::Instance();
        auto* rtdb = run->GetRuntimeDb();
        fNeulandHitPar = dynamic_cast<R3BNeulandHitPar*>(rtdb->findContainer(hitParName.c_str()));
        if (fNeulandHitPar != nullptr)
        {
            LOG(info) << "DigitizingTamex: HitPar " << hitParName
                      << " has been found in the root file. Using calibration values in rootfile.";
        }
        else
        {
            LOG(info) << "DigitizingTamex: HitPar " << hitParName << " cannot be found. Using default values.";
            fNeulandHitPar = nullptr;
        }
    }

    void Channel::AttachToPaddle(Digitizing::Paddle* /*v_paddle*/)
    {
        auto* paddle = GetPaddle();
        if (paddle == nullptr)
        {
            return;
        }

        if (CheckPaddleIDInHitPar())
        {
            fNeulandHitModulePar = fNeulandHitPar->GetModuleParAt(GetPaddle()->GetPaddleID() - 1);
            if (CheckPaddleIDInHitModulePar())
            {
                auto side = GetSide();
                par_.fSaturationCoefficient = fNeulandHitModulePar->GetPMTSaturation(static_cast<int>(side));
                par_.fEnergyGain = fNeulandHitModulePar->GetEnergyGain(static_cast<int>(side));
                par_.fPedestal = fNeulandHitModulePar->GetPedestal(static_cast<int>(side));
                par_.fPMTThresh = fNeulandHitModulePar->GetPMTThreshold(static_cast<int>(side));
                par_.fQdcMin = 1 / par_.fEnergyGain;
            }
        }
    }

    auto Channel::CheckPaddleIDInHitModulePar() const -> bool
    {
        auto is_valid = false;
        if (fNeulandHitModulePar == nullptr || GetPaddle() == nullptr)
        {
            return false;
        }

        if (GetPaddle()->GetPaddleID() != fNeulandHitModulePar->GetModuleId())
        {
            LOG(warn) << "Channel::SetHitModulePar:Wrong paddleID for the parameters!";
            is_valid = false;
        }
        else
        {
            is_valid = true;
        }
        return is_valid;
    }

    auto Channel::CheckPaddleIDInHitPar() const -> bool
    {
        auto is_valid = false;
        if (fNeulandHitPar == nullptr)
        {
            return false;
        }
        if (not fNeulandHitPar->hasChanged())
        {
            R3BLOG(warn, "Can't setup parameter in the root file correctly!.");
            return false;
        }

        auto PaddleId_max = fNeulandHitPar->GetNumModulePar();
        if (GetPaddle()->GetPaddleID() > PaddleId_max)
        {
            LOG(warn) << "Paddle id " << GetPaddle()->GetPaddleID() << " exceeds the id " << PaddleId_max
                      << " in the parameter file!";
            is_valid = false;
        }
        else
        {
            is_valid = true;
        }

        return is_valid;
    }

    void Channel::AddHit(Hit newHit)
    {
        if (newHit.time < par_.fTimeMin || newHit.time > par_.fTimeMax)
        {
            return;
        }
        InvalidateSignals();
        InvalidateTrigTime();
        fPMTPeaks.emplace_back(newHit, *this);
    }

    auto Channel::CreateSignal(const Peak& peak) const -> Signal
    {
        auto peakQdc = peak.GetQDC();
        auto peakTime = peak.GetLETime();
        auto qdc = ToQdc(peakQdc);

        auto signal = Signal{};
        signal.qdcUnSat = ToUnSatQdc(qdc);
        signal.qdc = qdc;
        signal.tdc = ToTdc(peakTime);
        signal.side = this->GetSide();
        LOG(debug) << "R3BDigitizingTamex: Create a signal with qdc " << signal.qdc << " and tdc " << signal.tdc
                   << std::endl;
        return signal;
    }

    template <typename Peak>
    void Channel::PeakPilingUp(/* inout */ std::vector<Peak>& peaks)
    {
        if (peaks.size() == 0)
        {
            return;
        }
        for (auto it = peaks.end() - 1; it != peaks.begin(); --it)
        {
            if (*it == *(it - 1))
            {
                *(it - 1) += *it;
                peaks.erase(it);
            }
        }
    }

    template <typename Peak>
    void Channel::ApplyThreshold(std::vector<Peak>& peaks)
    {
        // apply threshold on energy using c++ erase-remove idiom:
        auto it_end =
            std::remove_if(peaks.begin(),
                           peaks.end(),
                           [this](const auto& peak) { return peak.GetQDC() < this->GetParConstRef().fPMTThresh; });
        peaks.erase(it_end, peaks.end());
    }

    auto Channel::ConstructFQTPeaks(std::vector<PMTPeak>& pmtPeaks) -> std::vector<Peak>
    {
        auto FQTPeaks = std::vector<Peak>{};
        FQTPeaks.reserve(pmtPeaks.size());

        // sorting pmt peaks according to time:
        std::sort(pmtPeaks.begin(), pmtPeaks.end());

        PeakPilingUp(pmtPeaks);
        ApplyThreshold(pmtPeaks);
        for (auto const& peak : pmtPeaks)
        {
            FQTPeaks.emplace_back(peak, this);
        }
        return FQTPeaks;
    }

    auto Channel::ConstructSignals() -> Signals
    {
        fFQTPeaks = ConstructFQTPeaks(fPMTPeaks);
        // signal pileup:
        PeakPilingUp(fFQTPeaks);

        // construct Channel signals:
        auto signals = std::vector<Signal>{};
        signals.reserve(fFQTPeaks.size());

        for (const auto& peak : fFQTPeaks)
        {
            signals.emplace_back(CreateSignal(peak));
        }
        return signals;
    }

    auto Channel::GetFQTPeaks() -> const std::vector<Peak>&
    {

        if (!Is_ValidSignals())
        {
            ConstructSignals();
        }
        return fFQTPeaks;
    }

    auto Channel::GetPMTPeaks() -> const std::vector<PMTPeak>&
    {
        if (!Is_ValidSignals())
        {
            ConstructSignals();
        }
        return fPMTPeaks;
    }

    auto Channel::ToQdc(double qdc) const -> double
    {
        // apply energy smearing
        qdc = par_.fRnd->Gaus(qdc, par_.fEResRel * qdc);
        return qdc;
    }

    auto Channel::ToTdc(double time) const -> double { return time + par_.fRnd->Gaus(0., par_.fTimeRes); }

    auto Channel::ToUnSatQdc(double qdc) const -> double
    {
        // Apply reverse saturation
        if (par_.fExperimentalDataIsCorrectedForSaturation)
        {
            qdc = qdc / (1 - par_.fSaturationCoefficient * qdc);
        }
        // Apply reverse attenuation
        return qdc;
    }
} // namespace R3B::Digitizing::Neuland::Tamex
