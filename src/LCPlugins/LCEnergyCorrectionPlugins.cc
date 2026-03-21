/**
 *  @file   LCContent/src/LCPlugins/LCEnergyCorrectionPlugins.cc
 * 
 *  @brief  Implementation of the lc energy correction plugins class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCHelpers/ReclusterHelper.h"

#include "LCPlugins/LCEnergyCorrectionPlugins.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>

using namespace pandora;

namespace lc_content
{

namespace
{
using ClusterEnergyMap = std::unordered_map<const Cluster *, float>;

ClusterEnergyMap &GetSnapshotMap(const pandora::EnergyCorrectionType energyCorrectionType)
{
    static ClusterEnergyMap s_hadronicSnapshots;
    static ClusterEnergyMap s_electromagneticSnapshots;

    return (pandora::ELECTROMAGNETIC == energyCorrectionType) ? s_electromagneticSnapshots : s_hadronicSnapshots;
}
} // namespace

LCEnergyCorrectionPlugins::NonLinearityCorrection::NonLinearityCorrection(const FloatVector &inputEnergyCorrectionPoints,
        const FloatVector &outputEnergyCorrectionPoints) :
    m_inputEnergyCorrectionPoints(inputEnergyCorrectionPoints)
{
    const unsigned int nEnergyBins(m_inputEnergyCorrectionPoints.size());

    if (nEnergyBins != outputEnergyCorrectionPoints.size())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    for (unsigned int i = 0; i < nEnergyBins; ++i)
    {
        const float inputEnergy(m_inputEnergyCorrectionPoints.at(i));
        const float outputEnergy(outputEnergyCorrectionPoints.at(i));

        if (std::fabs(inputEnergy) < std::numeric_limits<float>::epsilon())
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

        m_energyCorrections.push_back(outputEnergy / inputEnergy);
    }

    if (nEnergyBins != m_energyCorrections.size())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode LCEnergyCorrectionPlugins::NonLinearityCorrection::MakeEnergyCorrections(const pandora::Cluster *const /*const pCluster*/, float &correctedEnergy) const
{
    const unsigned int nEnergyBins(m_energyCorrections.size());

    if (0 == nEnergyBins)
        return pandora::STATUS_CODE_SUCCESS;

    unsigned int index(nEnergyBins);

    for (unsigned int i = 0; i < nEnergyBins; ++i)
    {
        if (correctedEnergy < m_inputEnergyCorrectionPoints.at(i))
        {
            index = i;
            break;
        }
    }

    float correction(1.f);

    if ((0 == index) || (nEnergyBins == index))
    {
        correction = m_energyCorrections.at(std::min(index, nEnergyBins - 1));
    }
    else
    {
        const float lowCorrection(m_energyCorrections.at(index - 1)), highCorrection(m_energyCorrections.at(index));
        const float lowEnergy(m_inputEnergyCorrectionPoints.at(index - 1)), highEnergy(m_inputEnergyCorrectionPoints.at(index));
        correction = lowCorrection + (correctedEnergy - lowEnergy) * (highCorrection - lowCorrection) / (highEnergy - lowEnergy);
    }

    correctedEnergy *= correction;

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode LCEnergyCorrectionPlugins::NonLinearityCorrection::ReadSettings(const pandora::TiXmlHandle /*xmlHandle*/)
{
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LCEnergyCorrectionPlugins::CleanCluster::CleanCluster() :
    m_minCleanHitEnergy(0.5f),
    m_minCleanHitEnergyFraction(0.01f),
    m_minCleanCorrectedHitEnergy(0.1f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCEnergyCorrectionPlugins::CleanCluster::MakeEnergyCorrections(const Cluster *const pCluster, float &correctedHadronicEnergy) const
{
    const unsigned int firstPseudoLayer(this->GetPandora().GetPlugins()->GetPseudoLayerPlugin()->GetPseudoLayerAtIp());

    const float clusterHadronicEnergy(pCluster->GetHadronicEnergy());

    if (std::fabs(clusterHadronicEnergy) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    bool isFineGranularity(true);
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    // Loop over all constituent inner layer fine granularity hits, looking for anomalies
    for (OrderedCaloHitList::const_iterator layerIter = orderedCaloHitList.begin(), layerIterEnd = orderedCaloHitList.end();
        (layerIter != layerIterEnd) && isFineGranularity; ++layerIter)
    {
        const unsigned int pseudoLayer(layerIter->first);

        for (CaloHitList::const_iterator hitIter = layerIter->second->begin(), hitIterEnd = layerIter->second->end();
            hitIter != hitIterEnd; ++hitIter)
        {
            const CaloHit *const pCaloHit = *hitIter;

            if (this->GetPandora().GetGeometry()->GetHitTypeGranularity((*hitIter)->GetHitType()) > FINE)
            {
                isFineGranularity = false;
                break;
            }

            const float hitHadronicEnergy(pCaloHit->GetHadronicEnergy());

            if ((hitHadronicEnergy > m_minCleanHitEnergy) && (hitHadronicEnergy / clusterHadronicEnergy > m_minCleanHitEnergyFraction))
            {
                // Calculate new energy from surrounding layers
                float energyInPreviousLayer(0.);

                if (pseudoLayer > firstPseudoLayer)
                    energyInPreviousLayer = this->GetHadronicEnergyInLayer(orderedCaloHitList, pseudoLayer - 1);

                float energyInNextLayer(0.);

                if (pseudoLayer < std::numeric_limits<unsigned int>::max())
                    energyInNextLayer = this->GetHadronicEnergyInLayer(orderedCaloHitList, pseudoLayer + 1);

                const float energyInCurrentLayer = this->GetHadronicEnergyInLayer(orderedCaloHitList, pseudoLayer);

                // Calculate new energy estimate for hit and update cluster best energy estimate
                float energyInAdjacentLayers(energyInPreviousLayer + energyInNextLayer);

                if (pseudoLayer > firstPseudoLayer)
                    energyInAdjacentLayers /= 2.f;

                float newHitHadronicEnergy(energyInAdjacentLayers - energyInCurrentLayer + hitHadronicEnergy);
                newHitHadronicEnergy = std::max(newHitHadronicEnergy, m_minCleanCorrectedHitEnergy);

                if (newHitHadronicEnergy < hitHadronicEnergy)
                    correctedHadronicEnergy += newHitHadronicEnergy - hitHadronicEnergy;
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LCEnergyCorrectionPlugins::CleanCluster::GetHadronicEnergyInLayer(const OrderedCaloHitList &orderedCaloHitList, const unsigned int pseudoLayer) const
{
    OrderedCaloHitList::const_iterator iter = orderedCaloHitList.find(pseudoLayer);

    float hadronicEnergy(0.f);

    if (iter != orderedCaloHitList.end())
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            hadronicEnergy += (*hitIter)->GetHadronicEnergy();
        }
    }

    return hadronicEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCEnergyCorrectionPlugins::CleanCluster::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCleanHitEnergy", m_minCleanHitEnergy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCleanHitEnergyFraction", m_minCleanHitEnergyFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCleanCorrectedHitEnergy", m_minCleanCorrectedHitEnergy));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LCEnergyCorrectionPlugins::ScaleHotHadrons::ScaleHotHadrons() :
    m_minHitsForHotHadron(5),
    m_maxHitsForHotHadron(100),
    m_hotHadronInnerLayerCut(10),
    m_hotHadronMipFractionCut(0.4),
    m_hotHadronNHitsCut(50),
    m_hotHadronMipsPerHit(15.f),
    m_scaledHotHadronMipsPerHit(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCEnergyCorrectionPlugins::ScaleHotHadrons::MakeEnergyCorrections(const Cluster *const pCluster, float &correctedHadronicEnergy) const
{
    const unsigned int nHitsInCluster(pCluster->GetNCaloHits());

    // Initial hot hadron cuts
    if ((nHitsInCluster < m_minHitsForHotHadron) || (nHitsInCluster > m_maxHitsForHotHadron))
        return STATUS_CODE_SUCCESS;

    if ((pCluster->GetInnerPseudoLayer() < m_hotHadronInnerLayerCut) && (pCluster->GetMipFraction() < m_hotHadronMipFractionCut) &&
        (nHitsInCluster > m_hotHadronNHitsCut))
    {
        return STATUS_CODE_SUCCESS;
    }

    // Finally, check the number of mips per hit
    float clusterMipEnergy(0.);
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator layerIter = orderedCaloHitList.begin(), layerIterEnd = orderedCaloHitList.end();
        layerIter != layerIterEnd; ++layerIter)
    {
        for (CaloHitList::const_iterator hitIter = layerIter->second->begin(), hitIterEnd = layerIter->second->end();
            hitIter != hitIterEnd; ++hitIter)
        {
            clusterMipEnergy += (*hitIter)->GetMipEquivalentEnergy();
        }
    }

    const float meanMipsPerHit(clusterMipEnergy / static_cast<float>(nHitsInCluster));

    if ((meanMipsPerHit > 0.f) && (meanMipsPerHit > m_hotHadronMipsPerHit))
        correctedHadronicEnergy *= m_scaledHotHadronMipsPerHit / meanMipsPerHit;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCEnergyCorrectionPlugins::ScaleHotHadrons::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsForHotHadron", m_minHitsForHotHadron));

    if (0 == m_minHitsForHotHadron)
        return STATUS_CODE_INVALID_PARAMETER;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitsForHotHadron", m_maxHitsForHotHadron));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HotHadronInnerLayerCut", m_hotHadronInnerLayerCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HotHadronMipFractionCut", m_hotHadronMipFractionCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HotHadronNHitsCut", m_hotHadronNHitsCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HotHadronMipsPerHit", m_hotHadronMipsPerHit));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ScaledHotHadronMipsPerHit", m_scaledHotHadronMipsPerHit));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LCEnergyCorrectionPlugins::MuonCoilCorrection::MuonCoilCorrection() :
    m_muonHitEnergy(0.5f),
    m_coilEnergyLossCorrection(10.f),
    m_minMuonHitsInInnerLayer(3),
    m_coilEnergyCorrectionChi(3.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCEnergyCorrectionPlugins::MuonCoilCorrection::MakeEnergyCorrections(const Cluster *const pCluster, float &correctedHadronicEnergy) const
{
    bool containsMuonHit(false);
    unsigned int nMuonHitsInInnerLayer(0);
    unsigned int muonInnerLayer(std::numeric_limits<unsigned int>::max());

    // Extract muon-based properties from the cluster
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            if ((*hitIter)->GetHitType() == MUON)
            {
                containsMuonHit = true;
                ++nMuonHitsInInnerLayer;
            }
        }

        if (containsMuonHit)
        {
            muonInnerLayer = iter->first;
            break;
        }
    }

    if (!containsMuonHit)
        return STATUS_CODE_SUCCESS;;

    // Check whether energy deposits are likely to have been lost in coil region
    const CartesianVector muonInnerLayerCentroid(pCluster->GetCentroid(muonInnerLayer));
    const float centroidX(muonInnerLayerCentroid.GetX()), centroidY(muonInnerLayerCentroid.GetY());

    const float muonInnerLayerRadius(std::sqrt(centroidX * centroidX + centroidY * centroidY));
    const float coilInnerRadius(this->GetPandora().GetGeometry()->GetSubDetector(COIL).GetInnerRCoordinate());

    if (muonInnerLayerRadius < coilInnerRadius)
        return STATUS_CODE_SUCCESS;;

    const TrackList &trackList(pCluster->GetAssociatedTrackList());

    if (pCluster->GetInnerPseudoLayer() == muonInnerLayer)
    {
        // Energy correction for standalone muon cluster
        correctedHadronicEnergy += m_muonHitEnergy * static_cast<float>(nMuonHitsInInnerLayer);
    }
    else if (trackList.empty())
    {
        // Energy correction for neutral hadron cluster spilling into coil and muon detectors
        correctedHadronicEnergy += m_coilEnergyLossCorrection;
    }
    else
    {
        // Energy correction for charged hadron cluster spilling into coil and muon detectors
        if (nMuonHitsInInnerLayer < m_minMuonHitsInInnerLayer)
            return STATUS_CODE_SUCCESS;;

        float trackEnergySum(0.);

        for (TrackList::const_iterator iter = trackList.begin(), iterEnd = trackList.end(); iter != iterEnd; ++iter)
        {
            trackEnergySum += (*iter)->GetEnergyAtDca();
        }

        const float oldChi(ReclusterHelper::GetTrackClusterCompatibility(this->GetPandora(), correctedHadronicEnergy, trackEnergySum));
        const float newChi(ReclusterHelper::GetTrackClusterCompatibility(this->GetPandora(), correctedHadronicEnergy + m_coilEnergyLossCorrection, trackEnergySum));

        if ((oldChi < m_coilEnergyCorrectionChi) && (std::fabs(newChi) < std::fabs(oldChi)))
        {
            correctedHadronicEnergy += m_coilEnergyLossCorrection;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCEnergyCorrectionPlugins::MuonCoilCorrection::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MuonHitEnergy", m_muonHitEnergy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CoilEnergyLossCorrection", m_coilEnergyLossCorrection));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMuonHitsInInnerLayer", m_minMuonHitsInInnerLayer));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CoilEnergyCorrectionChi", m_coilEnergyCorrectionChi));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LCEnergyCorrectionPlugins::ThetaEnergyBinned::ThetaEnergyBinned(
    const FloatVector &ecalThetaBinEdges,
    const FloatVector &ecalEnergyBinEdges,
    const FloatVector &ecalScaleFactors,
    const FloatVector &hcalThetaBinEdges,
    const FloatVector &hcalEnergyBinEdges,
    const FloatVector &hcalScaleFactors,
    const pandora::EnergyCorrectionType energyCorrectionType) :
    m_energyCorrectionType(energyCorrectionType)
{
    m_ecalTable.m_thetaBinEdges = ecalThetaBinEdges;
    m_ecalTable.m_energyBinEdges = ecalEnergyBinEdges;
    m_ecalTable.m_scaleFactors = ecalScaleFactors;

    m_hcalTable.m_thetaBinEdges = hcalThetaBinEdges;
    m_hcalTable.m_energyBinEdges = hcalEnergyBinEdges;
    m_hcalTable.m_scaleFactors = hcalScaleFactors;

    ValidateTable(m_ecalTable);
    ValidateTable(m_hcalTable);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCEnergyCorrectionPlugins::ThetaEnergyBinned::MakeEnergyCorrections(const Cluster *const pCluster, float &correctedHadronicEnergy) const
{
    if (correctedHadronicEnergy <= std::numeric_limits<float>::epsilon())
        return STATUS_CODE_SUCCESS;

    ClusterEnergyMap &snapshotMap(GetSnapshotMap(m_energyCorrectionType));
    if (snapshotMap.end() == snapshotMap.find(pCluster))
        snapshotMap[pCluster] = correctedHadronicEnergy;

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.FillCaloHitList(caloHitList);
    caloHitList.insert(caloHitList.end(), pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());

    if (caloHitList.empty())
        return STATUS_CODE_SUCCESS;

    float ecalEnergy(0.f), hcalEnergy(0.f), totalWeight(0.f);
    CartesianVector weightedPosition(0.f, 0.f, 0.f);

    for (CaloHitList::const_iterator iter = caloHitList.begin(), iterEnd = caloHitList.end(); iter != iterEnd; ++iter)
    {
        const CaloHit *const pCaloHit(*iter);
        const float hitHadronicEnergy(std::max(0.f, pCaloHit->GetHadronicEnergy()));

        if (ECAL == pCaloHit->GetHitType())
            ecalEnergy += hitHadronicEnergy;
        else if (HCAL == pCaloHit->GetHitType())
            hcalEnergy += hitHadronicEnergy;

        weightedPosition += pCaloHit->GetPositionVector() * hitHadronicEnergy;
        totalWeight += hitHadronicEnergy;
    }

    if (totalWeight <= std::numeric_limits<float>::epsilon())
        return STATUS_CODE_SUCCESS;

    const CartesianVector clusterPosition(weightedPosition * (1.f / totalWeight));
    const float radius(clusterPosition.GetMagnitude());

    if (radius <= std::numeric_limits<float>::epsilon())
        return STATUS_CODE_SUCCESS;

    const float cosTheta(std::max(-1.f, std::min(1.f, clusterPosition.GetZ() / radius)));
    const float theta(std::acos(cosTheta));
    float correctedEcalEnergy(0.f), correctedHcalEnergy(0.f);

    if (ecalEnergy > std::numeric_limits<float>::epsilon())
        correctedEcalEnergy = ecalEnergy * LookupScale(m_ecalTable, theta, ecalEnergy);

    if (hcalEnergy > std::numeric_limits<float>::epsilon())
        correctedHcalEnergy = hcalEnergy * LookupScale(m_hcalTable, theta, hcalEnergy);

    correctedHadronicEnergy = correctedEcalEnergy + correctedHcalEnergy;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LCEnergyCorrectionPlugins::ThetaEnergyBinned::ResetLegacyEnergySnapshots()
{
    GetSnapshotMap(pandora::HADRONIC).clear();
    GetSnapshotMap(pandora::ELECTROMAGNETIC).clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LCEnergyCorrectionPlugins::ThetaEnergyBinned::GetLegacyEnergySnapshot(const Cluster *const pCluster,
    const pandora::EnergyCorrectionType energyCorrectionType, float &energy)
{
    const ClusterEnergyMap &snapshotMap(GetSnapshotMap(energyCorrectionType));
    const ClusterEnergyMap::const_iterator iter(snapshotMap.find(pCluster));

    if (snapshotMap.end() == iter)
        return false;

    energy = iter->second;
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCEnergyCorrectionPlugins::ThetaEnergyBinned::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LCEnergyCorrectionPlugins::ThetaEnergyBinned::FindBin(const FloatVector &edges, const float value)
{
    if (value < edges.front())
        return 0;

    if (value >= edges.back())
        return static_cast<unsigned int>(edges.size() - 2);

    const FloatVector::const_iterator upper(std::upper_bound(edges.begin(), edges.end(), value));
    return static_cast<unsigned int>(std::distance(edges.begin(), upper) - 1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LCEnergyCorrectionPlugins::ThetaEnergyBinned::LookupScale(const DomainTable &table, const float theta, const float energy)
{
    const unsigned int nEnergyBins(static_cast<unsigned int>(table.m_energyBinEdges.size() - 1));
    const unsigned int iTheta(FindBin(table.m_thetaBinEdges, theta));
    const unsigned int iEnergy(FindBin(table.m_energyBinEdges, energy));
    const unsigned int index(iTheta * nEnergyBins + iEnergy);

    return table.m_scaleFactors.at(index);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LCEnergyCorrectionPlugins::ThetaEnergyBinned::ValidateTable(const DomainTable &table)
{
    if ((table.m_thetaBinEdges.size() < 2) || (table.m_energyBinEdges.size() < 2))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if (!std::is_sorted(table.m_thetaBinEdges.begin(), table.m_thetaBinEdges.end()) ||
        !std::is_sorted(table.m_energyBinEdges.begin(), table.m_energyBinEdges.end()))
    {
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    for (std::size_t i = 1; i < table.m_thetaBinEdges.size(); ++i)
    {
        if (!(table.m_thetaBinEdges.at(i) > table.m_thetaBinEdges.at(i - 1)))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    for (std::size_t i = 1; i < table.m_energyBinEdges.size(); ++i)
    {
        if (!(table.m_energyBinEdges.at(i) > table.m_energyBinEdges.at(i - 1)))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    const std::size_t expectedSize((table.m_thetaBinEdges.size() - 1) * (table.m_energyBinEdges.size() - 1));
    if (table.m_scaleFactors.size() != expectedSize)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

} // namespace lc_content
