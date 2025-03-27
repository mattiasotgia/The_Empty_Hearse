#ifndef HELPER_H
#define HELPER_H

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/StandardRecord.h"

#include "TProfile.h"

#include <vector>     
#include <tuple>
#include <set>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>

const ana::SpillVar run    = ana::SIMPLESPILLVAR(hdr.run);
const ana::SpillVar subrun = ana::SIMPLESPILLVAR(hdr.subrun);
const ana::SpillVar event  = ana::SIMPLESPILLVAR(hdr.evt);

namespace logger {
    enum level {
        debug, 
        info, 
        warning, 
        error
    };

    std::map<level, std::string> log_level_string = {
        {debug, "debug"}, 
        {info, "info"}, 
        {warning, "warning"}, 
        {error, "error"}
    };

    std::ostream& log(level L, std::string where = "cout") {
        
        static std::ofstream file;

        std::ostream* writer = &std::cout;

        if (where != "cout") {
            file.open(where.c_str(), std::ios_base::app);
            writer = &file;
        }
        return *writer << "[" << log_level_string[L] << "] ";
    };
} // namespace log

namespace utils {
    int is_true_track(const caf::Proxy<caf::SRPFP> &pfp, const caf::SRSpillProxy *spill = nullptr) {
        if (pfp.trk.truth.p.pdg != pfp.shw.truth.p.pdg) {
            if (spill != nullptr)
                logger::log(logger::level::error) << "Found different pdg in run:event " << run(spill) << ":" << event(spill) << std::endl;
            else
                logger::log(logger::level::error) << "Found different pdg: " << pfp.trk.truth.p.pdg << " (pfp.trk.truth.p.pdg) and " << pfp.shw.truth.p.pdg << " (pfp.shw.truth.p.pdg)" << std::endl;
            return -1;
        }

        int true_pdg = (pfp.trk.truth.p.pdg == std::numeric_limits<int>::min()) ? pfp.shw.truth.p.pdg : pfp.trk.truth.p.pdg;

        std::vector<int> track_pdgs{
            13, 211, 2212
        };

        return std::find(track_pdgs.begin(), track_pdgs.end(), std::abs(true_pdg)) != track_pdgs.end();
    };  
    
    template<class T>
    const bool sanity_check(T object) {

        // Check normal types
        if (std::is_same_v<T, int> || std::is_same_v<T, double> || std::is_same_v<T, float>) {
            if (std::isnan(object)) return false;
        }

        // Check vectors
        if (std::is_same_v<T, caf::SRVector3D>) {
            if (std::isnan(object.x) || std::isnan(object.y) || std::isnan(object.z)) return false;
        }

        return true;
    }
} // namespace utils


namespace bins {
    const ana::Binning simple = ana::Binning::Simple(15, 0, 1);
    const ana::Binning simple_log = ana::Binning::LogUniform(50, 0, 1);
} // namespace bins

namespace hard_code {

    // std::vector<unsigned> bad_events {
    //     10201, 10202, 10203, 10204, 10205, 10206, 10207, 10208, 10209, 10210, 10211, 10212, 10213, 10214, 10215, 10216, 10217, 10218, 10219, 10220, 10221, 10222, 10223, 10224, 10225, 10226, 10227, 10228, 10229, 10230, 10231, 10232, 10233, 10234, 10235, 10236, 10237, 10238, 10239, 10240, 10241, 10242, 10243, 10244, 10245, 10246, 10247, 10248, 10249, 10250, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 12551, 12552, 12553, 12554, 12555, 12556, 12557, 12558, 12559, 12560, 12561, 12562, 12563, 12564, 12565, 12566, 12567, 12568, 12569, 12570, 12571, 12572, 12573, 12574, 12575, 12576, 12577, 12578, 12579, 12580, 12581, 12582, 12583, 12584, 12585, 12586, 12587, 12588, 12589, 12590, 12591, 12592, 12593, 12594, 12595, 12596, 12597, 12598, 12599, 12600, 14951, 14952, 14953, 14954, 14955, 14956, 14957, 14958, 14959, 14960, 14961, 14962, 14963, 14964, 14965, 14966, 14967, 14968, 14969, 14970, 14971, 14972, 14973, 14974, 14975, 14976, 14977, 14978, 14979, 14980, 14981, 14982, 14983, 14984, 14985, 14986, 14987, 14988, 14989, 14990, 14991, 14992, 14993, 14994, 14995, 14996, 14997, 14998, 14999, 15000
    // };

    std::vector<unsigned> bad_events = {
        10151, 10152, 10153, 10154, 10155, 10156, 10157, 10158, 10159, 10160, 10161, 10162, 10163, 10164, 10165, 10166, 10167, 10168, 10169, 10170, 10171, 10172, 10173, 10174, 10175, 10176, 10177, 10178, 10179, 10180, 10181, 10182, 10183, 10184, 10185, 10186, 10187, 10188, 10189, 10190, 10191, 10192, 10193, 10194, 10195, 10196, 10197, 10198, 10199, 10200, 3951, 3952, 3953, 3954, 3955, 3956, 3957, 3958, 3959, 3960, 3961, 3962, 3963, 3964, 3965, 3966, 3967, 3968, 3969, 3970, 3971, 3972, 3973, 3974, 3975, 3976, 3977, 3978, 3979, 3980, 3981, 3982, 3983, 3984, 3985, 3986, 3987, 3988, 3989, 3990, 3991, 3992, 3993, 3994, 3995, 3996, 3997, 3998, 3999, 4000, 14801, 14802, 14803, 14804, 14805, 14806, 14807, 14808, 14809, 14810, 14811, 14812, 14813, 14814, 14815, 14816, 14817, 14818, 14819, 14820, 14821, 14822, 14823, 14824, 14825, 14826, 14827, 14828, 14829, 14830, 14831, 14832, 14833, 14834, 14835, 14836, 14837, 14838, 14839, 14840, 14841, 14842, 14843, 14844, 14845, 14846, 14847, 14848, 14849, 14850, 10701, 10702, 10703, 10704, 10705, 10706, 10707, 10708, 10709, 10710, 10711, 10712, 10713, 10714, 10715, 10716, 10717, 10718, 10719, 10720, 10721, 10722, 10723, 10724, 10725, 10726, 10727, 10728, 10729, 10730, 10731, 10732, 10733, 10734, 10735, 10736, 10737, 10738, 10739, 10740, 10741, 10742, 10743, 10744, 10745, 10746, 10747, 10748, 10749, 10750, 15201, 15202, 15203, 15204, 15205, 15206, 15207, 15208, 15209, 15210, 15211, 15212, 15213, 15214, 15215, 15216, 15217, 15218, 15219, 15220, 15221, 15222, 15223, 15224, 15225, 15226, 15227, 15228, 15229, 15230, 15231, 15232, 15233, 15234, 15235, 15236, 15237, 15238, 15239, 15240, 15241, 15242, 15243, 15244, 15245, 15246, 15247, 15248, 15249, 15250, 15901, 15902, 15903, 15904, 15905, 15906, 15907, 15908, 15909, 15910, 15911, 15912, 15913, 15914, 15915, 15916, 15917, 15918, 15919, 15920, 15921, 15922, 15923, 15924, 15925, 15926, 15927, 15928, 15929, 15930, 15931, 15932, 15933, 15934, 15935, 15936, 15937, 15938, 15939, 15940, 15941, 15942, 15943, 15944, 15945, 15946, 15947, 15948, 15949, 15950, 13651, 13652, 13653, 13654, 13655, 13656, 13657, 13658, 13659, 13660, 13661, 13662, 13663, 13664, 13665, 13666, 13667, 13668, 13669, 13670, 13671, 13672, 13673, 13674, 13675, 13676, 13677, 13678, 13679, 13680, 13681, 13682, 13683, 13684, 13685, 13686, 13687, 13688, 13689, 13690, 13691, 13692, 13693, 13694, 13695, 13696, 13697, 13698, 13699, 13700, 1901, 1902, 1903, 1904, 1905, 1906, 1907, 1908, 1909, 1910, 1911, 1912, 1913, 1914, 1915, 1916, 1917, 1918, 1919, 1920, 1921, 1922, 1923, 1924, 1925, 1926, 1927, 1928, 1929, 1930, 1931, 1932, 1933, 1934, 1935, 1936, 1937, 1938, 1939, 1940, 1941, 1942, 1943, 1944, 1945, 1946, 1947, 1948, 1949, 1950, 11101, 11102, 11103, 11104, 11105, 11106, 11107, 11108, 11109, 11110, 11111, 11112, 11113, 11114, 11115, 11116, 11117, 11118, 11119, 11120, 11121, 11122, 11123, 11124, 11125, 11126, 11127, 11128, 11129, 11130, 11131, 11132, 11133, 11134, 11135, 11136, 11137, 11138, 11139, 11140, 11141, 11142, 11143, 11144, 11145, 11146, 11147, 11148, 11149, 11150, 5701, 5702, 5703, 5704, 5705, 5706, 5707, 5708, 5709, 5710, 5711, 5712, 5713, 5714, 5715, 5716, 5717, 5718, 5719, 5720, 5721, 5722, 5723, 5724, 5725, 5726, 5727, 5728, 5729, 5730, 5731, 5732, 5733, 5734, 5735, 5736, 5737, 5738, 5739, 5740, 5741, 5742, 5743, 5744, 5745, 5746, 5747, 5748, 5749, 5750, 8601, 8602, 8603, 8604, 8605, 8606, 8607, 8608, 8609, 8610, 8611, 8612, 8613, 8614, 8615, 8616, 8617, 8618, 8619, 8620, 8621, 8622, 8623, 8624, 8625, 8626, 8627, 8628, 8629, 8630, 8631, 8632, 8633, 8634, 8635, 8636, 8637, 8638, 8639, 8640, 8641, 8642, 8643, 8644, 8645, 8646, 8647, 8648, 8649, 8650, 1251, 1252, 1253, 1254, 1255, 1256, 1257, 1258, 1259, 1260, 1261, 1262, 1263, 1264, 1265, 1266, 1267, 1268, 1269, 1270, 1271, 1272, 1273, 1274, 1275, 1276, 1277, 1278, 1279, 1280, 1281, 1282, 1283, 1284, 1285, 1286, 1287, 1288, 1289, 1290, 1291, 1292, 1293, 1294, 1295, 1296, 1297, 1298, 1299, 1300, 10551, 10552, 10553, 10554, 10555, 10556, 10557, 10558, 10559, 10560, 10561, 10562, 10563, 10564, 10565, 10566, 10567, 10568, 10569, 10570, 10571, 10572, 10573, 10574, 10575, 10576, 10577, 10578, 10579, 10580, 10581, 10582, 10583, 10584, 10585, 10586, 10587, 10588, 10589, 10590, 10591, 10592, 10593, 10594, 10595, 10596, 10597, 10598, 10599, 10600, 10001, 10002, 10003, 10004, 10005, 10006, 10007, 10008, 10009, 10010, 10011, 10012, 10013, 10014, 10015, 10016, 10017, 10018, 10019, 10020, 10021, 10022, 10023, 10024, 10025, 10026, 10027, 10028, 10029, 10030, 10031, 10032, 10033, 10034, 10035, 10036, 10037, 10038, 10039, 10040, 10041, 10042, 10043, 10044, 10045, 10046, 10047, 10048, 10049, 10050, 15601, 15602, 15603, 15604, 15605, 15606, 15607, 15608, 15609, 15610, 15611, 15612, 15613, 15614, 15615, 15616, 15617, 15618, 15619, 15620, 15621, 15622, 15623, 15624, 15625, 15626, 15627, 15628, 15629, 15630, 15631, 15632, 15633, 15634, 15635, 15636, 15637, 15638, 15639, 15640, 15641, 15642, 15643, 15644, 15645, 15646, 15647, 15648, 15649, 15650, 5901, 5902, 5903, 5904, 5905, 5906, 5907, 5908, 5909, 5910, 5911, 5912, 5913, 5914, 5915, 5916, 5917, 5918, 5919, 5920, 5921, 5922, 5923, 5924, 5925, 5926, 5927, 5928, 5929, 5930, 5931, 5932, 5933, 5934, 5935, 5936, 5937, 5938, 5939, 5940, 5941, 5942, 5943, 5944, 5945, 5946, 5947, 5948, 5949, 5950, 9201, 9202, 9203, 9204, 9205, 9206, 9207, 9208, 9209, 9210, 9211, 9212, 9213, 9214, 9215, 9216, 9217, 9218, 9219, 9220, 9221, 9222, 9223, 9224, 9225, 9226, 9227, 9228, 9229, 9230, 9231, 9232, 9233, 9234, 9235, 9236, 9237, 9238, 9239, 9240, 9241, 9242, 9243, 9244, 9245, 9246, 9247, 9248, 9249, 9250, 14501, 14502, 14503, 14504, 14505, 14506, 14507, 14508, 14509, 14510, 14511, 14512, 14513, 14514, 14515, 14516, 14517, 14518, 14519, 14520, 14521, 14522, 14523, 14524, 14525, 14526, 14527, 14528, 14529, 14530, 14531, 14532, 14533, 14534, 14535, 14536, 14537, 14538, 14539, 14540, 14541, 14542, 14543, 14544, 14545, 14546, 14547, 14548, 14549, 14550, 17301, 17302, 17303, 17304, 17305, 17306, 17307, 17308, 17309, 17310, 17311, 17312, 17313, 17314, 17315, 17316, 17317, 17318, 17319, 17320, 17321, 17322, 17323, 17324, 17325, 17326, 17327, 17328, 17329, 17330, 17331, 17332, 17333, 17334, 17335, 17336, 17337, 17338, 17339, 17340, 17341, 17342, 17343, 17344, 17345, 17346, 17347, 17348, 17349, 17350, 7901, 7902, 7903, 7904, 7905, 7906, 7907, 7908, 7909, 7910, 7911, 7912, 7913, 7914, 7915, 7916, 7917, 7918, 7919, 7920, 7921, 7922, 7923, 7924, 7925, 7926, 7927, 7928, 7929, 7930, 7931, 7932, 7933, 7934, 7935, 7936, 7937, 7938, 7939, 7940, 7941, 7942, 7943, 7944, 7945, 7946, 7947, 7948, 7949, 7950, 16101, 16102, 16103, 16104, 16105, 16106, 16107, 16108, 16109, 16110, 16111, 16112, 16113, 16114, 16115, 16116, 16117, 16118, 16119, 16120, 16121, 16122, 16123, 16124, 16125, 16126, 16127, 16128, 16129, 16130, 16131, 16132, 16133, 16134, 16135, 16136, 16137, 16138, 16139, 16140, 16141, 16142, 16143, 16144, 16145, 16146, 16147, 16148, 16149, 16150, 11451, 11452, 11453, 11454, 11455, 11456, 11457, 11458, 11459, 11460, 11461, 11462, 11463, 11464, 11465, 11466, 11467, 11468, 11469, 11470, 11471, 11472, 11473, 11474, 11475, 11476, 11477, 11478, 11479, 11480, 11481, 11482, 11483, 11484, 11485, 11486, 11487, 11488, 11489, 11490, 11491, 11492, 11493, 11494, 11495, 11496, 11497, 11498, 11499, 11500, 16401, 16402, 16403, 16404, 16405, 16406, 16407, 16408, 16409, 16410, 16411, 16412, 16413, 16414, 16415, 16416, 16417, 16418, 16419, 16420, 16421, 16422, 16423, 16424, 16425, 16426, 16427, 16428, 16429, 16430, 16431, 16432, 16433, 16434, 16435, 16436, 16437, 16438, 16439, 16440, 16441, 16442, 16443, 16444, 16445, 16446, 16447, 16448, 16449, 16450, 10801, 10802, 10803, 10804, 10805, 10806, 10807, 10808, 10809, 10810, 10811, 10812, 10813, 10814, 10815, 10816, 10817, 10818, 10819, 10820, 10821, 10822, 10823, 10824, 10825, 10826, 10827, 10828, 10829, 10830, 10831, 10832, 10833, 10834, 10835, 10836, 10837, 10838, 10839, 10840, 10841, 10842, 10843, 10844, 10845, 10846, 10847, 10848, 10849, 10850, 11001, 11002, 11003, 11004, 11005, 11006, 11007, 11008, 11009, 11010, 11011, 11012, 11013, 11014, 11015, 11016, 11017, 11018, 11019, 11020, 11021, 11022, 11023, 11024, 11025, 11026, 11027, 11028, 11029, 11030, 11031, 11032, 11033, 11034, 11035, 11036, 11037, 11038, 11039, 11040, 11041, 11042, 11043, 11044, 11045, 11046, 11047, 11048, 11049, 11050, 14351, 14352, 14353, 14354, 14355, 14356, 14357, 14358, 14359, 14360, 14361, 14362, 14363, 14364, 14365, 14366, 14367, 14368, 14369, 14370, 14371, 14372, 14373, 14374, 14375, 14376, 14377, 14378, 14379, 14380, 14381, 14382, 14383, 14384, 14385, 14386, 14387, 14388, 14389, 14390, 14391, 14392, 14393, 14394, 14395, 14396, 14397, 14398, 14399, 14400, 2351, 2352, 2353, 2354, 2355, 2356, 2357, 2358, 2359, 2360, 2361, 2362, 2363, 2364, 2365, 2366, 2367, 2368, 2369, 2370, 2371, 2372, 2373, 2374, 2375, 2376, 2377, 2378, 2379, 2380, 2381, 2382, 2383, 2384, 2385, 2386, 2387, 2388, 2389, 2390, 2391, 2392, 2393, 2394, 2395, 2396, 2397, 2398, 2399, 2400, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 2751, 2752, 2753, 2754, 2755, 2756, 2757, 2758, 2759, 2760, 2761, 2762, 2763, 2764, 2765, 2766, 2767, 2768, 2769, 2770, 2771, 2772, 2773, 2774, 2775, 2776, 2777, 2778, 2779, 2780, 2781, 2782, 2783, 2784, 2785, 2786, 2787, 2788, 2789, 2790, 2791, 2792, 2793, 2794, 2795, 2796, 2797, 2798, 2799, 2800, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749, 750, 8751, 8752, 8753, 8754, 8755, 8756, 8757, 8758, 8759, 8760, 8761, 8762, 8763, 8764, 8765, 8766, 8767, 8768, 8769, 8770, 8771, 8772, 8773, 8774, 8775, 8776, 8777, 8778, 8779, 8780, 8781, 8782, 8783, 8784, 8785, 8786, 8787, 8788, 8789, 8790, 8791, 8792, 8793, 8794, 8795, 8796, 8797, 8798, 8799, 8800, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 9051, 9052, 9053, 9054, 9055, 9056, 9057, 9058, 9059, 9060, 9061, 9062, 9063, 9064, 9065, 9066, 9067, 9068, 9069, 9070, 9071, 9072, 9073, 9074, 9075, 9076, 9077, 9078, 9079, 9080, 9081, 9082, 9083, 9084, 9085, 9086, 9087, 9088, 9089, 9090, 9091, 9092, 9093, 9094, 9095, 9096, 9097, 9098, 9099, 9100, 14201, 14202, 14203, 14204, 14205, 14206, 14207, 14208, 14209, 14210, 14211, 14212, 14213, 14214, 14215, 14216, 14217, 14218, 14219, 14220, 14221, 14222, 14223, 14224, 14225, 14226, 14227, 14228, 14229, 14230, 14231, 14232, 14233, 14234, 14235, 14236, 14237, 14238, 14239, 14240, 14241, 14242, 14243, 14244, 14245, 14246, 14247, 14248, 14249, 14250, 8051, 8052, 8053, 8054, 8055, 8056, 8057, 8058, 8059, 8060, 8061, 8062, 8063, 8064, 8065, 8066, 8067, 8068, 8069, 8070, 8071, 8072, 8073, 8074, 8075, 8076, 8077, 8078, 8079, 8080, 8081, 8082, 8083, 8084, 8085, 8086, 8087, 8088, 8089, 8090, 8091, 8092, 8093, 8094, 8095, 8096, 8097, 8098, 8099, 8100, 10101, 10102, 10103, 10104, 10105, 10106, 10107, 10108, 10109, 10110, 10111, 10112, 10113, 10114, 10115, 10116, 10117, 10118, 10119, 10120, 10121, 10122, 10123, 10124, 10125, 10126, 10127, 10128, 10129, 10130, 10131, 10132, 10133, 10134, 10135, 10136, 10137, 10138, 10139, 10140, 10141, 10142, 10143, 10144, 10145, 10146, 10147, 10148, 10149, 10150, 12551, 12552, 12553, 12554, 12555, 12556, 12557, 12558, 12559, 12560, 12561, 12562, 12563, 12564, 12565, 12566, 12567, 12568, 12569, 12570, 12571, 12572, 12573, 12574, 12575, 12576, 12577, 12578, 12579, 12580, 12581, 12582, 12583, 12584, 12585, 12586, 12587, 12588, 12589, 12590, 12591, 12592, 12593, 12594, 12595, 12596, 12597, 12598, 12599, 12600, 12701, 12702, 12703, 12704, 12705, 12706, 12707, 12708, 12709, 12710, 12711, 12712, 12713, 12714, 12715, 12716, 12717, 12718, 12719, 12720, 12721, 12722, 12723, 12724, 12725, 12726, 12727, 12728, 12729, 12730, 12731, 12732, 12733, 12734, 12735, 12736, 12737, 12738, 12739, 12740, 12741, 12742, 12743, 12744, 12745, 12746, 12747, 12748, 12749, 12750, 10201, 10202, 10203, 10204, 10205, 10206, 10207, 10208, 10209, 10210, 10211, 10212, 10213, 10214, 10215, 10216, 10217, 10218, 10219, 10220, 10221, 10222, 10223, 10224, 10225, 10226, 10227, 10228, 10229, 10230, 10231, 10232, 10233, 10234, 10235, 10236, 10237, 10238, 10239, 10240, 10241, 10242, 10243, 10244, 10245, 10246, 10247, 10248, 10249, 10250, 13551, 13552, 13553, 13554, 13555, 13556, 13557, 13558, 13559, 13560, 13561, 13562, 13563, 13564, 13565, 13566, 13567, 13568, 13569, 13570, 13571, 13572, 13573, 13574, 13575, 13576, 13577, 13578, 13579, 13580, 13581, 13582, 13583, 13584, 13585, 13586, 13587, 13588, 13589, 13590, 13591, 13592, 13593, 13594, 13595, 13596, 13597, 13598, 13599, 13600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650
    };

    const bool is_bad_event(const caf::SRSpillProxy *sr) {
        return std::find(bad_events.begin(), bad_events.end(), static_cast<int>(event(sr))) != bad_events.end();
    };
} // namespace hard_code

namespace cheating {

    struct chi2
    {
        double muon;
        double proton;
        double kaon;
        double pi;
    }; 

    struct values_minmax {
        double min;
        double max;
    }; 

    double dist_cut = 10.;
    values_minmax barycenterFM_deltaZ_Trigger = {0., 100.};

    std::string dEdx_temp = 
        "/exp/icarus/app/users/msotgia/analysis/sbnana_v09_93_01_thesis_analysis/analysis/dEdxrestemplates.root";
    TFile* file = TFile::Open(dEdx_temp.c_str());

    auto dedx_range_pro = (TProfile *)file->Get("dedx_range_pro");
    auto dedx_range_ka  = (TProfile *)file->Get("dedx_range_ka");
    auto dedx_range_pi  = (TProfile *)file->Get("dedx_range_pi");
    auto dedx_range_mu  = (TProfile *)file->Get("dedx_range_mu");

    chi2 chi2_ALG (std::vector<double> &dEdx, std::vector<double> &RR, double rr_min, double rr_max) {
        /* The output is chi2s
         * Those are used to classify final state particles 
        */

        double threshold = 0.5;
        double max_rr = rr_max; // Max value for the residual range (RR)
        double min_rr = rr_min; // Min value for the residual range (RR)

        std::vector<float> trkdedx;
        std::vector<float> trkres;
        std::vector<double> vpida;

        for (std::size_t i(0); i < dEdx.size(); ++i) {
            if (i == 0 || i == dEdx.size() - 1)
                continue;
            if (RR[i] < max_rr && RR[i] > rr_min) {
                trkdedx.push_back(dEdx[i]);
                trkres.push_back(RR[i]);
            }
        }

        int npt = 0;
        double chi2pro = 0;
        double chi2ka = 0;
        double chi2pi = 0;
        double chi2mu = 0;
        double avgdedx = 0;
        double PIDA = 0;

        int used_trkres = 0;
        for (unsigned i = 0; i < trkdedx.size(); ++i) { 
            // hits
            // ignore the first and the last point
            // if (i == 0 || i == trkdedx.size() - 1) continue;
            avgdedx += trkdedx[i];
            if (trkres[i] < 26) {
                PIDA += trkdedx[i] * std::pow(trkres[i], 0.42);
                vpida.push_back(trkdedx[i] * std::pow(trkres[i], 0.42));
                used_trkres++;
            }
            if (trkdedx[i] > 100 || trkdedx[i] < threshold)
                continue; // protect against large pulse height
            
            int bin = dedx_range_pro->FindBin(trkres[i]);
            if (bin >= 1 && bin <= dedx_range_pro->GetNbinsX()) {
                
                double bincpro = dedx_range_pro->GetBinContent(bin);
                if (bincpro < 1e-6)  // for 0 bin content, using neighboring bins
                    bincpro = 
                        (dedx_range_pro->GetBinContent(bin - 1) + dedx_range_pro->GetBinContent(bin + 1)) / 2;
                
                double bincka = dedx_range_ka->GetBinContent(bin);
                if (bincka < 1e-6)
                    bincka =
                        (dedx_range_ka->GetBinContent(bin - 1) + dedx_range_ka->GetBinContent(bin + 1)) / 2;

                double bincpi = dedx_range_pi->GetBinContent(bin);
                if (bincpi < 1e-6)
                    bincpi =
                        (dedx_range_pi->GetBinContent(bin - 1) + dedx_range_pi->GetBinContent(bin + 1)) / 2;
                
                double bincmu = dedx_range_mu->GetBinContent(bin);
                if (bincmu < 1e-6)
                    bincmu =
                        (dedx_range_mu->GetBinContent(bin - 1) + dedx_range_mu->GetBinContent(bin + 1)) / 2;
                
                double binepro = dedx_range_pro->GetBinError(bin);
                if (binepro < 1e-6)
                    binepro =
                        (dedx_range_pro->GetBinError(bin - 1) + dedx_range_pro->GetBinError(bin + 1)) / 2;

                double bineka = dedx_range_ka->GetBinError(bin);
                if (bineka < 1e-6)
                    bineka = (dedx_range_ka->GetBinError(bin - 1) + dedx_range_ka->GetBinError(bin + 1)) / 2;

                double binepi = dedx_range_pi->GetBinError(bin);
                if (binepi < 1e-6)
                    binepi = (dedx_range_pi->GetBinError(bin - 1) + dedx_range_pi->GetBinError(bin + 1)) / 2;

                double binemu = dedx_range_mu->GetBinError(bin);
                if (binemu < 1e-6)
                    binemu = (dedx_range_mu->GetBinError(bin - 1) + dedx_range_mu->GetBinError(bin + 1)) / 2;

                // double errke = 0.05*trkdedx[i];   //5% KE resolution

                double errdedx = 0.04231 + 0.0001783 * trkdedx[i] * trkdedx[i]; // resolution on dE/dx
                
                errdedx *= trkdedx[i];
                
                chi2pro += std::pow((trkdedx[i] - bincpro) / std::sqrt(std::pow(binepro, 2) + std::pow(errdedx, 2)), 2);
                chi2ka  += std::pow((trkdedx[i] - bincka)  / std::sqrt(std::pow(bineka, 2)  + std::pow(errdedx, 2)), 2);
                chi2pi  += std::pow((trkdedx[i] - bincpi)  / std::sqrt(std::pow(binepi, 2)  + std::pow(errdedx, 2)), 2);
                chi2mu  += std::pow((trkdedx[i] - bincmu)  / std::sqrt(std::pow(binemu, 2)  + std::pow(errdedx, 2)), 2);
                // std::cout<<i<<" "<<trkdedx[i]<<" "<<trkres[i]<<" "<<bincpro<<std::endl;
                ++npt;
            }
        } // hits

        return {chi2mu / npt, chi2pro / npt, chi2ka / npt, chi2pi / npt};
    } // chi2 chi2_ALG

    using level_t = logger::level;
 
    // // 
    // // Definition of a chi2 statistics to test the cheat level (of the reconstruction goodness)
    // template<class T>
    // const T chi2(T __a, T __b, std::string label = "undef") {
    //     double chi2 = 0; // default "as if truthmatched";
    //     if (__b == 0 ) {
    //         if (__a == __b) {
    //             logger::log(level_t::warning) << "Found truth value zero (checking label " 
    //                                           << label << "), but reco == truth, so this will be ignored"
    //                                           << std::endl;
    //             return chi2;
    //         } else {
    //             logger::log(level_t::warning) << "Found truth value zero (checking label " 
    //                                           << label << "), but reco != truth, so returning their difference squared..."
    //                                           << std::endl;
    //             return std::pow(std::abs(__a - __b), 2);
    //         }
    //     }
    //     chi2 = std::pow(std::abs(__a - __b)/(double)__b, 2);
    //     logger::log(level_t::info) << label << " chi2: " << chi2 << "(__a = " << __a << "; __b = " << __b << ")" << std::endl;
    //     return chi2;
    // };


    // // The idea is to build a strong estimate (chi2-like) to test the goodness of 
    // // this cheating... This is integrated over all the variables
    // // At the moment the only possibility is ( V_reco - V_truth ) / V_truth
    // // Minimal selection applied: iscc + pdg == 14 + event is not in the bad_list...
    // const ana::MultiVar score ([](const caf::SRSliceProxy *slice) -> std::vector<double> {

    //     std::vector<double> score;

    //     score.push_back(chi2<double>(slice->vertex.x, slice->truth.position.x, "vertex_position"));
    //     score.push_back(chi2<double>(slice->vertex.y, slice->truth.position.y, "vertex_position"));
    //     score.push_back(chi2<double>(slice->vertex.z, slice->truth.position.z, "vertex_position"));

    //     // recob::PFParticle and other stuff
    //     for (auto const &pfp: slice->reco.pfp) {
    //         if (utils::is_true_track(pfp) == -1) continue;

    //         bool is_track = utils::is_true_track(pfp);

    //         // Track score as additional score value...
    //         // score.push_back(chi2<double>(pfp.trackScore, is_track ? 1 : 0, "track_score"));
            
    //         if (is_track) {
    //             score.push_back(chi2<double>(pfp.trk.truth.bestmatch.hit_completeness, 1, "hit_completeness"));
    //             score.push_back(chi2<double>(pfp.trk.truth.bestmatch.hit_purity, 1, "hit_purity"));

    //             score.push_back(chi2<double>(pfp.trk.start.x, pfp.trk.truth.p.start.x, "track_start_point"));
    //             score.push_back(chi2<double>(pfp.trk.start.y, pfp.trk.truth.p.start.y, "track_start_point"));
    //             score.push_back(chi2<double>(pfp.trk.start.z, pfp.trk.truth.p.start.z, "track_start_point"));

    //             score.push_back(chi2<double>(pfp.trk.end.x, pfp.trk.truth.p.end.x, "track_end_point"));
    //             score.push_back(chi2<double>(pfp.trk.end.y, pfp.trk.truth.p.end.y, "track_end_point"));
    //             score.push_back(chi2<double>(pfp.trk.end.z, pfp.trk.truth.p.end.z, "track_end_point")); 

    //             score.push_back(chi2<double>(pfp.trk.len, pfp.trk.truth.p.length, "track_lenght"));

    //             score.push_back(chi2<double>(
    //                 slice->vertex.x-pfp.trk.start.x, 
    //                 slice->truth.position.x-pfp.trk.truth.p.start.x, "vertex_start_position"
    //             ));
    //             score.push_back(chi2<double>(
    //                 slice->vertex.y-pfp.trk.start.y, 
    //                 slice->truth.position.y-pfp.trk.truth.p.start.y, "vertex_start_position"
    //             ));
    //             score.push_back(chi2<double>(
    //                 slice->vertex.z-pfp.trk.start.z, 
    //                 slice->truth.position.z-pfp.trk.truth.p.start.z, "vertex_start_position"
    //             ));
    //         } else {
    //             score.push_back(chi2<double>(pfp.shw.truth.bestmatch.hit_completeness, 1));
    //             score.push_back(chi2<double>(pfp.shw.truth.bestmatch.hit_purity, 1));

    //             score.push_back(chi2<double>(pfp.shw.start.x, pfp.shw.truth.p.start.x, "shower_start_point"));
    //             score.push_back(chi2<double>(pfp.shw.start.y, pfp.shw.truth.p.start.y, "shower_start_point"));
    //             score.push_back(chi2<double>(pfp.shw.start.z, pfp.shw.truth.p.start.z, "shower_start_point"));

    //             score.push_back(chi2<double>(pfp.shw.end.x, pfp.shw.truth.p.end.x, "shower_end_point"));
    //             score.push_back(chi2<double>(pfp.shw.end.y, pfp.shw.truth.p.end.y, "shower_end_point"));
    //             score.push_back(chi2<double>(pfp.shw.end.z, pfp.shw.truth.p.end.z, "shower_end_point")); 
                
    //             score.push_back(chi2<double>(pfp.shw.len, pfp.shw.truth.p.length, "shower_lenght"));

    //             score.push_back(chi2<double>(
    //                 slice->vertex.x-pfp.shw.start.x,
    //                 slice->truth.position.x-pfp.shw.truth.p.start.x, "vertex_start_position"
    //             ));
    //             score.push_back(chi2<double>(
    //                 slice->vertex.y-pfp.shw.start.y,
    //                 slice->truth.position.y-pfp.shw.truth.p.start.y, "vertex_start_position"
    //             ));
    //             score.push_back(chi2<double>(
    //                 slice->vertex.z-pfp.shw.start.z,
    //                 slice->truth.position.z-pfp.shw.truth.p.start.z, "vertex_start_position"
    //             ));
    //         }

    //     } // loop pfp
    //     return score;
    // }); // score ana::MultiVar

    std::map<bool, std::string> boolean_print = {
        {true, "true"},
        {false, "false"}
    };

    const ana::SpillMultiVar test_variables ([](const caf::SRSpillProxy *spill) -> std::vector<double> {
        std::vector<double> tested_variables;

        std::cout << "This is a spill, run:event = " << run(spill) << ":" << event(spill) << std::endl;

        std::cout << "Looking at truth interaction..." << std::endl;
        std::cout << "--> Found nnu = " << spill->mc.nnu << " true nu(s)" << std::endl;

        int inu = 0;
        for (auto const& nu: spill->mc.nu) {

            std::cout << "Nu no. " << inu << ": pdg = "  << nu.pdg << ", initpdg = " << nu.initpdg << ", index = " << nu.index << std::endl;
            
            std::cout << "          isnc = " << boolean_print[nu.isnc] << ", "
                      << "iscc = " << boolean_print[nu.iscc] << ", "
                      << "isvtxcont = " << boolean_print[nu.isvtxcont] << ", "
                      << "is_numucc_primary = " << boolean_print[nu.is_numucc_primary] << ", "
                      <<std::endl; 
            
            // std::cout << "          vtx = ("
            //           << nu.vtx.x << ", " 
            //           << nu.vtx.y << ", " 
            //           << nu.vtx.z << "), " << std::endl;

            std::cout << "          position = ("
                      << nu.position.x << ", " 
                      << nu.position.y << ", " 
                      << nu.position.z << "), " << std::endl;

            std::cout << "          npiplus = " << nu.npiplus << ", "
                      << "npiminus = " << nu.npiminus << ", "
                      << "npizero = " << nu.npizero << ", "
                      << "nproton = " << nu.nproton << ", "
                      << "nneutron = " << nu.nneutron << " (BEFORE final state interaction!)"
                      << std::endl;
            
            std::cout << "This nu has nprim = " << nu.nprim << " primary particle(s) " << std::endl;
            std::cout << "The primaries pf this interaction have pdg in {" << std::endl;

            for (auto const& prim: nu.prim) 
                std::cout << "\t" << prim.pdg << ", \tlenght = " << prim.length << " cm" << std::endl;

            std::cout << "}" << std::endl;

            inu ++ ;
        } // loop over nu(s)


        std::cout << "Looking at slices..." << std::endl;
        std::cout << "--> Found nslc = " << spill->nslc << " slice(s) " << std::endl;
        
        int islc = 0;
        for (auto const& slice: spill->slc) {

            std::cout << "Slice no. " << islc << " with ID = " << slice.self << std::endl;
            std::cout << "This is a " << ((slice.truth.iscc && slice.truth.pdg == 14) ? "GOOD" : "BAD") << " slice" << std::endl;

            std::cout << "Reco vertex in (" 
                      << slice.vertex.x << ", " 
                      << slice.vertex.y << ", " 
                      << slice.vertex.z << ") and true vertex in ("
                      << slice.truth.position.x << ", " 
                      << slice.truth.position.y << ", " 
                      << slice.truth.position.z << ")" 
                      << std::endl;
            std::cout << "This slice has npfp = " << slice.reco.npfp << " pfp(s)" << std::endl;
            std::cout << "Looping though primaries = {" << std::flush;

            for (auto const& iprimary: slice.primary) std::cout << iprimary << ", " << std::flush;

            std::cout << "} " << std::endl;

            int ipfp = 0;
            for (auto const& pfp: slice.reco.pfp) {

                // compute new chi2
                std::vector<double> dedx;
                std::vector<double> rr;
                int use_plane = 2;
                for (std::size_t ihit = 0; ihit < pfp.trk.calo[use_plane].points.size(); ++ihit) {
                    dedx.push_back(pfp.trk.calo[use_plane].points[ihit].dedx);
                    rr.push_back(pfp.trk.calo[use_plane].points[ihit].rr);
                } // calo points

                // input to chi2_ALG vector dedx, vector rr, rr_min, rr_max
                // output chi2s {chi2mu/npt,chi2pro/npt,chi2ka/npt,chi2pi/npt}
                chi2 chi2_values = chi2_ALG(dedx, rr, 0.0, 25.0);
                
                std::cout << "[ipfp = " << ipfp << "] --> " 
                          << "ID = " << pfp.id << ", "
                          << "parent = " << pfp.parent << ", "
                          << "parent_is_primary = " << boolean_print[pfp.parent_is_primary] << ", "
                          << "trackScore = " << pfp.trackScore << ", "
                          << "(for checking purpouses, .pfochar.vtxdist = " << pfp.pfochar.vtxdist << "), "
                          << " chi2 (proton hypotesis) = " << chi2_values.proton << ", "
                          << " chi2 (muon hypotesis) = " << chi2_values.muon << ", "
                          << "slcID = " << pfp.slcID << std::endl;
                
                if (pfp.parent == -1) continue; // Neutrino PFP :) all other variables are empty

                if (utils::is_true_track(pfp, spill) == true) {

                    ///////////////////////////////////////////
                    ////////////////// TRACK //////////////////
                    ///////////////////////////////////////////

                                //[ipfp = 0]
                    std::cout << "           --> TRK properties: " << std::endl;
                    std::cout << "                             .trk:         " 
                              << "len = " << pfp.trk.len << " cm, "
                              << "start = (" << pfp.trk.start.x << ", " << pfp.trk.start.y << ", " <<pfp.trk.start.z << "), "
                              << "end = (" << pfp.trk.end.x << ", " << pfp.trk.end.y << ", " <<pfp.trk.end.z << "), "
                              << std::endl
                              << "                                           "
                              << "dir = (" << pfp.trk.dir.x << ", " << pfp.trk.dir.y << ", " <<pfp.trk.dir.z << "), "
                              << "dir_end = (" << pfp.trk.dir_end.x << ", " << pfp.trk.dir_end.y << ", " <<pfp.trk.dir_end.z << "), "
                              << "costh = " << pfp.trk.costh << ", "
                              << "phi = " << pfp.trk.phi << ", "
                              << std::endl;

                    std::cout << "                             .trk.truth.p: "
                              << "len = " << pfp.trk.truth.p.length << " cm, "
                              << "start = (" << pfp.trk.truth.p.start.x << ", " << pfp.trk.truth.p.start.y << ", " <<pfp.trk.truth.p.start.z << "), "
                              << "end = (" << pfp.trk.truth.p.end.x << ", " << pfp.trk.truth.p.end.y << ", " <<pfp.trk.truth.p.end.z << "), "
                              << "gen = (" << pfp.trk.truth.p.gen.x << ", " << pfp.trk.truth.p.gen.y << ", " <<pfp.trk.truth.p.gen.z << "), "
                              << std::endl
                              << "                                           "
                              << "pdg = " << pfp.trk.truth.p.pdg << ", "
                              << "cont_tpc = " << boolean_print[pfp.trk.truth.p.cont_tpc] << ", "
                              << "crosses_tpc = " << boolean_print[pfp.trk.truth.p.crosses_tpc] << ", "
                              << "contained = " << boolean_print[pfp.trk.truth.p.contained] << ", "
                              << "cryostat = " << pfp.trk.truth.p.cryostat
                              << std::endl;
                    std::cout << "                             .trk.truth.bestmatch: "
                              << "hit_completeness = " << pfp.trk.truth.bestmatch.hit_completeness << ", "
                              << "hit_purity = " << pfp.trk.truth.bestmatch.hit_purity
                              << std::endl;
                } else {

                    ///////////////////////////////////////////
                    ///////////////// SHOWER //////////////////
                    ///////////////////////////////////////////
                    
                    std::cout << "           --> SHW properties: " << std::endl;
                    std::cout << "                             .shw:         " 
                              << "len = " << pfp.shw.len << " cm, "
                              << "start = (" << pfp.shw.start.x << ", " << pfp.shw.start.y << ", " <<pfp.shw.start.z << "), "
                              << "end = (" << pfp.shw.end.x << ", " << pfp.shw.end.y << ", " <<pfp.shw.end.z << "), "
                              << std::endl
                              << "                                           "
                              << "dir = (" << pfp.shw.dir.x << ", " << pfp.shw.dir.y << ", " <<pfp.shw.dir.z << "), "
                              << "conversion_gap = " << pfp.shw.conversion_gap << " cm, "
                              << "density = " << pfp.shw.density << " MeV/cm, "
                              << "open_angle = " << pfp.shw.open_angle << " rad"
                              << std::endl;
                    std::cout << "                             .shw.truth.p: "
                              << "len = " << pfp.shw.truth.p.length << " cm, "
                              << "start = (" << pfp.shw.truth.p.start.x << ", " << pfp.shw.truth.p.start.y << ", " <<pfp.shw.truth.p.start.z << "), "
                              << "end = (" << pfp.shw.truth.p.end.x << ", " << pfp.shw.truth.p.end.y << ", " <<pfp.shw.truth.p.end.z << "), "
                              << "gen = (" << pfp.shw.truth.p.gen.x << ", " << pfp.shw.truth.p.gen.y << ", " <<pfp.shw.truth.p.gen.z << "), "
                              << std::endl
                              << "                                           "
                              << "pdg = " << pfp.shw.truth.p.pdg << ", "
                              << "cont_tpc = " << boolean_print[pfp.shw.truth.p.cont_tpc] << ", "
                              << "crosses_tpc = " << boolean_print[pfp.shw.truth.p.crosses_tpc] << ", "
                              << "contained = " << boolean_print[pfp.shw.truth.p.contained] << ", "
                              << "cryostat = " << pfp.shw.truth.p.cryostat
                              << std::endl;
                    std::cout << "                             .shw.truth.bestmatch: "
                              << "hit_completeness = " << pfp.shw.truth.bestmatch.hit_completeness << ", "
                              << "hit_purity = " << pfp.shw.truth.bestmatch.hit_purity
                              << std::endl;
                }

                ipfp ++ ;
            } // loop over pfp(s)

            islc ++ ;
        } // loop over slice(s)

        return tested_variables;
    });

    const ana::Cut cut_numuCC ([](const caf::SRSliceProxy *slice) -> bool {
        logger::log(level_t::info) << "This is a " << ((slice->truth.iscc && slice->truth.pdg == 14) ? "GOOD" : "BAD") << " slice" << std::endl;
        return slice->truth.iscc && slice->truth.pdg == 14;
    }); 

    const ana::SpillCut cut_bad_events([](const caf::SRSpillProxy *spill) -> bool {
        // logger::log(level_t::info) << "This is a " << (hard_code::is_bad_event(spill) ? "BAD" : "GOOD") << " run:event" << std::endl;
        return !hard_code::is_bad_event(spill);
    });

} // namespace cheating

#endif // HELPER_H