#include "google/s2/s1angle.h"
#include "google/s2/s2cellid.h"
#include "google/s2/s2latlng.h"
#include "google/s2/s2cellunion.h"
#include "google/s2/s2regioncoverer.h"
#include "google/s2/s2cap.h"
#include "radar.h"
#include "matrix.h"
#include <unordered_map>
#include <map>

static const double kEarthRadius = 6378137;
static const double kEarthCircumference = 40075017;
static const int    kGroundLevel = 23;

inline double Distance(const DPoint &a, const DPoint &b) {
    return kEarthRadius * a.id.ToLatLng().GetDistance(b.id.ToLatLng()).radians();
}

inline double Distance(const DPoint &a, S2LatLng &b) {
    return kEarthRadius * a.id.ToLatLng().GetDistance(b).radians();
}

// deprecated, use GetCloseLevel from s2
inline int ClosetLevel(double d) {
    return d < 1 ? kGroundLevel : kGroundLevel - std::min(kGroundLevel-1, (int)ceil(log2(d)));
}

inline double EarthMetersToRadians(double meters) {
    return (2 * M_PI) * (meters/kEarthCircumference); 
}

inline double RadiansToEarthMeters(double radians) {
    return (radians * kEarthCircumference) / (2 * M_PI);
}

// Returns the cell level with a side most closely matching the
// specified number of meters.
static int GetClosestLevel(double meters) {
    return S2::kAvgEdge.GetClosestLevel(EarthMetersToRadians(meters));
}

inline S1Angle ClosesetAngle(double dc) {
    int level = GetClosestLevel(dc);
    S2CellId sid = S2CellId::FromLatLng(S2LatLng::FromDegrees(39,116)).parent(level);
    S2Cap cap = S2Cell(sid).GetCapBound();
    return cap.angle();
}

// 用cannopy的思想对空间点集快速聚类
int quick_union(const std::vector<DPoint> &points, std::vector<Cluster> &cls, double dc) {
    int level = GetClosestLevel(dc);
    // cellid -> #cells, cluster
    std::unordered_map<S2CellId, std::pair<int, int> > grids;
    std::unordered_map<S2CellId, std::pair<int, int> >::iterator iter;

    for(int i = 0, l = points.size(); i < l; ++ i) {
        S2CellId big_grid_id = points[i].id.parent(level); 
        iter = grids.find(big_grid_id);
        if (iter == grids.end()) 
            grids[big_grid_id] = std::make_pair(points[i].samples, -1);
        else
            grids[big_grid_id].first += points[i].samples;
    }

    int nr_center = 0;
    std::vector<S2CellId> neibors;
    std::unordered_map<S2CellId, std::pair<int,int> >::iterator nb_iter;
    for(iter = grids.begin(); iter != grids.end(); ++iter) {
        if (iter->second.second == -1) {
            iter->second.second = nr_center;
            nr_center ++;
        } else continue; 

        neibors.clear();
        iter->first.AppendAllNeighbors(level, &neibors);
        for(size_t i = 0; i < neibors.size(); ++ i) {
            S2CellId nb = neibors[i];
            nb_iter = grids.find(nb);
            if (nb_iter != grids.end() && nb_iter->second.second == -1) {
                nb_iter->second.second = iter->second.second;
                neibors.push_back(nb_iter->first);
                nb_iter->first.AppendAllNeighbors(level, &neibors);
            }
        }
        // printf("%d->%d\n", nr_center, neibors.size());
    }

    cls.resize(nr_center);
    for(iter = grids.begin(); iter != grids.end(); ++iter) {
        int cnt = iter->second.first;
        int clz = iter->second.second;
        S2LatLng ll = iter->first.ToLatLng();
        cls[clz].lat += cnt * ll.lat().degrees();
        cls[clz].lng += cnt * ll.lng().degrees();
        cls[clz].cnt += cnt;
    }
    for(auto &clz: cls) {
        if (clz.cnt == 0) {
            clz.lat = 0;
            clz.lng = 0;
        } else {
            clz.lat /= clz.cnt;
            clz.lng /= clz.cnt;
        }
    }
    return cls.size();
}


// see: http://en.wikipedia.org/wiki/DBSCAN 
int dbscan(std::vector<DPoint> &points, double dc, int minPts) {
    int level = GetClosestLevel(dc);
    typedef std::unordered_map<S2CellId, std::vector<int> > hashmap;
    hashmap cell_index;
    for(int i = 0, l = points.size(); i < l; ++ i) {
        S2CellId pid = points[i].id.parent(level);
        hashmap::iterator it = cell_index.find(pid);
        if (it == cell_index.end()) {
            std::vector<int> cells;
            cells.push_back(i);
            cell_index.insert(std::make_pair(pid, cells));
        } else it->second.push_back(i);
    }

    std::vector<S2CellId> neibors;
    std::vector<int> candidates;
    std::vector<int> pending;

    neibors.reserve(8);
    candidates.reserve(points.size());
    pending.reserve(points.size());

    bool visited[points.size()];
    bool visiting[points.size()];
    for(size_t i = 0; i < points.size(); ++ i) {
        visited[i] = visiting[i] = false;
    }
    // memset(&visited[0], points.size()*sizeof(bool), 0);
    // memset(&visiting[0], points.size()*sizeof(bool), 0);

    int clazz = 1, main_cluster = 0, main_cluster_size = 0;
    for(int c = 0, l = points.size(); c < l; ++ c) {
        if (visited[c]) continue;
        // walk the graph from points[c]
        pending.clear();
        pending.push_back(c);

        for(size_t k = 0; k < pending.size(); ++k) {
            int i = pending[k];
            if (visited[i]) continue;
            visited[i] = true;
            // printf("#%d[%lu]:", i, pending.size());

            neibors.clear();
            candidates.clear();

            S2CellId pid = points[i].id.parent(level);
            neibors.push_back(pid);
            pid.AppendAllNeighbors(pid.level(), &neibors);
            for(S2CellId nb: neibors) {
                hashmap::iterator it = cell_index.find(nb);
                if (it != cell_index.end()) {
                    for(int idx: it->second) {
                        candidates.push_back(idx);
                    }
                }
            }

            // 邻节点可能被多次加入,需要消重
            std::sort(candidates.begin(), candidates.end());
            std::vector<int>::iterator end_it;
            end_it = std::unique(candidates.begin(), candidates.end());
            candidates.erase(end_it, candidates.end());

            // 按最大距离过滤出真正的邻居节点
            S2LatLng l1 = points[i].id.ToLatLng();
            for(std::vector<int>::iterator cand_it = candidates.begin(); 
                    cand_it != candidates.end(); ) {
                S2LatLng l2 = points[*cand_it].id.ToLatLng();
                if (kEarthRadius * l1.GetDistance(l2).radians() > dc) {
                    cand_it = candidates.erase(cand_it);
                } else {
                    cand_it ++;
                }
            }
            // 划分核心节点和边缘节点
            points[i].is_core = false;
            if (points[i].samples - 1 > minPts || (int)candidates.size() > minPts) {
                points[i].is_core = true;
            } else {
                int real_nb_cnt = 0;
                for(auto &c: candidates) {
                    real_nb_cnt += points[c].samples;
                }
                if (real_nb_cnt > minPts) {
                    points[i].is_core = true;
                }
            }
            
            if (points[i].is_core) {
                for(int j: candidates) {
                    // printf(" %d-%d", j, visiting[j]);
                    // avoid add Vj to pending list twice
                    if (!visiting[j]) {
                        pending.push_back(j);
                        visiting[j] = true;
                    }
                }
            }
        }
        
        int pending_samples = 0;
        for(int i: pending) {
            points[i].cluster = clazz;
            pending_samples += points[i].samples;
        }
        if (main_cluster_size < (int)pending.size()) {
            main_cluster_size = pending_samples;
            main_cluster = clazz;
        }
        ++ clazz;
    }

    int cluster_size_counter[clazz];
    for(int i = 0; i < clazz; ++ i) {
        cluster_size_counter[i] = 0;
    }
    for(const auto &pt: points) {
        cluster_size_counter[pt.cluster] += pt.samples;
    }
    std::vector<std::pair<int,int>> remark;
    for(int i = 0; i < clazz; ++ i) {
        remark.push_back(std::make_pair(cluster_size_counter[i], i));
    }
    std::sort(remark.begin(), remark.end(),
            [](const std::pair<int, int>&a, const std::pair<int, int>&b){return a.first > b.first; }
            );
    for(int c = 0; c < clazz; ++ c) {
        cluster_size_counter[remark[c].second] = c;
    }
    
    // 俺类大小重新分配标签, 0 为最大的类
    for(auto &pt: points) {
        pt.cluster = cluster_size_counter[pt.cluster];
    }

    return clazz;
}

// wrap function for dbscan
int dbscan(std::vector<S2CellId> cells, std::vector<int>& cluster, double dc, int minPts) {
    std::vector<DPoint> points;
    for(auto cid: cells) {
        S2LatLng ll = cid.ToLatLng();
        double lat = ll.lat().degrees();
        double lng = ll.lat().degrees();
        points.push_back(DPoint(lat, lng));
    }
    int clazz = dbscan(points, dc, minPts);
    cluster.resize(clazz);
    for(int i = 0, l = cells.size(); i < l; ++ i) {
        cluster.push_back(points[i].cluster);
    }
    return clazz;
    
}

int dbscan(double *lat, double *lng, short *cluster, int n, double dc, int minPts) {
    std::vector<DPoint> points;
    for(int i = 0; i < n; ++ i) {
        points.push_back(DPoint(lat[i], lng[i]));
    }
    dbscan(points, dc, minPts);
    for(int i = 0; i < n; ++ i) {
        cluster[i] = (points[i].cluster << 1) + (int)points[i].is_core;
    }
    return n;
}

// 对点集进行快速聚类,O(n)的复杂度
int density_peak_cluster_fast(std::vector<DPoint>& points, double dc) {
    if (points.empty()) return 0; 
    if (points.size() == 1) {
        points[0].cluster = 0;
        return 1;
    }
    std::sort(points.begin(), points.end());
    points[0].cluster = 0;
    int clazz = 0, density = 0, best_clazz = 0, max_density = 1;
    for(int i = 1, l = points.size(); i < l; ++ i) {
        int j = i - 1;
        DPoint &pi = points[i];
        DPoint &pj = points[j];
        float d = Distance(pi, pj);
        if (d > dc) {
            if (density > max_density) {
                best_clazz = clazz;
                max_density = density;
            }
            ++ clazz;
            density = 1;
        }
        pi.cluster = clazz; 
    }
    // 最后一个点可能跟第一个点相近
    if (points[0].cluster != clazz && dc > Distance(points.back(), points[0])) {
        for(int i = points.size() - 1; i >= 0 && points[i].cluster == clazz; --i) {
            points[i].cluster = points[0].cluster;
        }
        clazz --;
    }
    return clazz;
}

// 快速计算密度最大的点
S2LatLng density_peak_centroid(std::vector<DPoint>& points, double dc) {
    int max_clazz = 1 + density_peak_cluster_fast(points, dc);
    int count[max_clazz];
    memset(&count[0], 0, sizeof(count[0]) * max_clazz);
    for(int i = 0, l = points.size(); i < l; ++ i) {
        count[points[i].cluster] += 1;
        // printf("%d -> C-%d\n", i, points[i].clazz);
    }
    int max_count = 0, best_clazz = 0;
    for(int i = 0; i < max_clazz; ++ i) {
        if (count[i] > max_count) {
            max_count = count[i];
            best_clazz = i;
        }
    }
    int total = 0;
    double lat = 0, lng = 0;
    for(int i = 0, l = points.size(); i < l; ++ i) {
        if (points[i].cluster == best_clazz) {
            lat += points[i].lat;
            lng += points[i].lng;
            ++ total;
        }
    }
    if (total > 0) {
        lat /= total;
        lng /= total;
    }
    return S2LatLng::FromDegrees(lat, lng);
}

// 快速过滤密度低的区域
int density_peak_filter(DPointVec& points, double threshold, double dc) {
    S2LatLng cent = density_peak_centroid(points, dc);
    // printf("cent = %.6f,%.6f, thres=%g/%g\n", lat, lng, threshold, dc);
    typedef std::vector<DPoint>::iterator dpit;
    for(dpit it = points.begin(); it != points.end(); ++it) {
        double d = Distance(*it, cent);
        it->is_noise = d > threshold;
    }
    return points.size();
}

//把level设为1m
void ground_to_1m(DPointVec& points) {
    std::vector<DPoint>::iterator pit;
    for(pit = points.begin(); pit != points.end(); ++ pit) {
        pit->id = pit->id.parent(kGroundLevel);
    }
}

// 点密度估计. 点密度 = 距离<dc的邻居样本数 + 自身样本数
void construct_density_map(DPointVec& points,
        double dc, std::unordered_map<S2CellId, int> &density_map) {
    std::vector<S2CellId> neibors;
    std::unordered_map<S2CellId, int> counter;
    std::unordered_map<S2CellId, int>::iterator it;
    for(const auto &pt: points) {
        // ignore the noise
        if (pt.is_noise)
            continue;
        S2CellId cid(pt.id);
        it = counter.find(cid);
        if (it == counter.end()) {
            counter.insert(std::make_pair(cid, pt.samples));
        } else {
            it->second += pt.samples;
        }
    }

    S2RegionCoverer cover;
    cover.set_max_level(kGroundLevel);
    cover.set_min_level(kGroundLevel);
    S1Angle angle = ClosesetAngle(dc);

    for(it = counter.begin(); it != counter.end(); ++it) {
        S2Cap cap = S2Cap::FromAxisAngle(it->first.ToPoint(), angle) ;
        neibors.clear();
        cover.GetCovering(cap, &neibors);

        int count = it->second;
        for(auto c: neibors) {
            auto cnt_iter = counter.find(c);
            if (cnt_iter != counter.end()) {
                count += cnt_iter->second;
            }
            // std::cout<<*nb_it<<std::endl;
        }
        density_map[it->first] = count;
    } 
}

void fill_density(DPointVec& points, std::unordered_map<S2CellId, int> density_map) {
    for(auto &pt: points) {
        S2CellId sid(pt.id);
        auto it = density_map.find(sid);
        if (it != density_map.end()) 
            pt.density = it->second;
        else
            pt.density = 0;
    } 
}

// 用最大密度估计点集的中心
void estimate_density(DPointVec& points, double dc) {
    if (points.empty()) return;
    // 缩放到1m大小的网格
    if (points[0].id.level() > kGroundLevel)
        ground_to_1m(points);
    std::unordered_map<S2CellId, int> density_map;
    construct_density_map(points, dc, density_map);
    // 保存计算结果 
    fill_density(points, density_map);
}

int density_peak_find_centroid(DPointVec& points, double dc) {
    estimate_density(points, dc);
    int max_density = 0, centroid = 0;
    for(int i = 0, l = points.size(); i < l; ++ i) {
        if(points[i].density > max_density) {
            max_density = points[i].density;
            centroid = i;
        }
    }
    return centroid;
}

S2LatLng density_peak_centroid_slow(DPointVec& points, double dc) {
    int centroid = density_peak_find_centroid(points, dc);
    return points[centroid].id.ToLatLng();
}

int density_peak_filter_slow(DPointVec& points, double threshold, double dc) {
    S2LatLng cent = density_peak_centroid_slow(points, dc);
    std::vector<DPoint>::iterator pit;
    for(pit = points.begin(); pit != points.end(); ++ pit) {
        double d = Distance(*pit, cent);
        pit->is_noise = d > threshold;
    }
    // printf("centroid: %.6f %.6f %g/%g c=%d d=%d n=%d\n", lat, lng, threshold, centroid, max_density, points.size());
    return points.size();
}

int dpfilt(double *lat, double *lng, int n, double threshold, double dc) {
    std::vector<DPoint> points;
    for(int i = 0; i < n; ++ i) {
        points.push_back(DPoint(lat[i], lng[i]));
    }
    // printf("n = %d, t=%g/%g\n", n, threshold, dc);
    n = density_peak_filter(points, threshold, dc);
    for(int i = 0; i < n; ++ i) {
        lat[i] = points[i].lat;
        lng[i] = points[i].lng;
    }
    return n;
}

int dpfilt_slow(double *lat, double *lng, int n, double threshold, double dc) {
    std::vector<DPoint> points;
    for(int i = 0; i < n; ++ i) {
        points.push_back(DPoint(lat[i], lng[i]));
    }
    // printf("n = %d, t=%g/%g\n", n, threshold, dc);
    n = density_peak_filter_slow(points, threshold, dc);
    for(int i = 0; i < n; ++ i) {
        lat[i] = points[i].lat;
        lng[i] = points[i].lng;
    }
    return n;
}

bool has_multi_center(DPointVec &points, double min_expand_radius, double max_expand_radius) {
    std::vector<Cluster> unions; 
    int nr_center = 0;

    if (points.size() < 50) {
        nr_center = quick_union(points, unions, max_expand_radius);
        if (nr_center > 1)
            return true;
        unions.clear();
    }
    
    nr_center = quick_union(points, unions, min_expand_radius);
    int nr_robust_center = 0;
    int limit = round(points.size()*0.03);
    for(auto &c: unions) {
        if (c.cnt > limit)
            nr_robust_center ++;
    }
    // printf("center = %d, points = %d\n", nr_robust_center, points.size());
    return nr_robust_center > 1;
}

void mark_noise(DPointVec &points, double ratio) {
    int samples = 0;
    for(const auto &pt: points) {
        samples += pt.samples;
    }
    int min_cluster_size = samples * ratio + 0.5;
    if (min_cluster_size == 0) return ;
    std::map<int, int> labels;
    for (auto pt: points) {
        if (labels.find(pt.cluster) == labels.end()) {
            labels[pt.cluster] = pt.samples;
        } else {
            labels[pt.cluster] += pt.samples;
        }
    }
    for (auto &pt: points) {
        if(labels[pt.cluster] < min_cluster_size) {
            pt.is_noise = true;
        }
    }
}

int scan_noise(DPointVec &points,
        double dbscan_min_pts,
        double dbscan_min_radius,
        double noise_level) {
    int nr_cluster = dbscan(points, dbscan_min_radius, dbscan_min_pts);
    // keep top3 cluster
    if (nr_cluster > 3) {
        for(auto &pt: points) {
            if (pt.cluster > 2) {
                pt.is_noise = true;
            }
        }
    }

    mark_noise(points, noise_level);
    int nr_noise = 0;
    for(auto &pt: points) {
        if (pt.is_noise) nr_noise ++;
    }
    return nr_noise;
}

// 删除所有的噪音点
void prune(DPointVec &points, PruneOption opt) {
    scan_noise(points, opt.dbscan_min_radius, opt.dbscan_min_pts, opt.noise_level); 
    std::vector<DPoint> output;
    std::copy_if(points.begin(), points.end(), std::back_insert_iterator<DPointVec>(output),
            [](const DPoint &p1) { return p1.is_noise; });
    points.swap(output);
}

// 判断是否是移动wifi
bool is_mobile(DPointVec &points, MobileDetectOption opt) {
    // printf("==>%d\n", points.size());
    if (points.size() < 2) return false;
    if (has_multi_center(points, opt.min_expand_radius, opt.max_expand_radius))
        return true;

    if (points.size() >= (size_t)opt.min_pts) {
        int nr_points = points.size();
        int nr_noise  = scan_noise(points, opt.dbscan_min_pts, opt.dbscan_min_radius, opt.noise_level);
        int main_cluster_size = 0;
        for(int i = 0, l = points.size(); i < l; ++ i) {
            if (points[i].cluster == 0) {
                main_cluster_size ++;
            } 
        }

        if (nr_noise > main_cluster_size || main_cluster_size*opt.core_scale < nr_points) 
            return true;
    }
    return false;
}

bool maybe_intersect(
        const std::unordered_set<S2CellId> &c0,
        DPointVec &points, int level, int c) {
    std::vector<S2CellId> neibors;
    for(const auto &pt: points) {
        S2CellId pid = pt.id.parent(level);
        if (pt.cluster == c) {
            if (c0.find(pid) != c0.end())
                return true;
            pid.AppendAllNeighbors(level, &neibors);
            for(auto nid: neibors) {
                if (c0.find(nid) != c0.end())
                    return true;
            }
            neibors.clear();
        }
    }
    return false;
}

void merge_cluster_0(DPointVec &points, double distance) {
    int level = GetClosestLevel(distance);
    std::unordered_map<S2CellId, int> c0;
    for(const auto &pt: points) {
        if (pt.cluster == 0) {
            S2CellId pid = pt.id.parent(level);
            auto iter = c0.find(pid);
            if (iter == c0.end()) {
                // c0.insert(std::make_pair(pid, pt.samples));
                c0.insert(std::make_pair(pid, 1));
            } else {
                // iter->second += pt.samples;
                iter->second ++;
            }
        }
    }

    std::vector<S2CellId> neibors;
    for(auto &pt: points) {
        if (pt.cluster > 2)
            continue;
        S2CellId pid = pt.id.parent(level);
        neibors.push_back(pid);
        pid.AppendAllNeighbors(level, &neibors);
        int nr_neighbors = 0;
        for(auto nid: neibors) {
            auto iter = c0.find(nid);
            if (iter != c0.end()) {
                nr_neighbors += iter->second;
            } 
        }
        neibors.clear();
        if (nr_neighbors > 1) {
            pt.cluster = 0;
        }
    }
}

void merge_cluster_1(DPointVec &points, double distance) {
    int level = GetClosestLevel(distance);
    std::unordered_set<S2CellId> c0;
    int cluster_size[6] = {0, 0, 0, 0, 0, 0};
    for(const auto &pt: points) {
        if (pt.cluster < 5) {
            cluster_size[pt.cluster] ++;
            if (pt.cluster == 0) 
                c0.insert(pt.id.parent(level));
        }
    }
    
    for(int c = 1; c < 3; ++ c) {
        if (cluster_size[c] < 6)
            break;
        if(maybe_intersect(c0, points, level, c)) {
            for(auto &pt: points) {
                if (pt.cluster == c) {
                    pt.cluster = 0;
                }
            }
        }
    }
}

static inline bool float_equal(double x, double y) {
    return fabs(x-y) < 1e-16;
}

bool compare_dpoint(const DPoint &d1, const DPoint &d2) {
    if (!float_equal(d1.lng, d2.lng)) {
        return d1.lat < d2.lat;
    }
    return d1.lng < d2.lng;
}

struct DPointAngleComparator {
    DPointAngleComparator(double x, double y):
        orig_x(x),
        orig_y(y)
    { }

    bool operator()(const DPoint &d1, const DPoint &d2) {
        double dx1 = d1.lng - orig_x;
        double dy1 = d1.lat - orig_y;
        double dx2 = d2.lng - orig_x;
        double dy2 = d2.lat - orig_y;

        double norm1 = sqrt(dx1*dx1 + dy1*dy1);
        double norm2 = sqrt(dx2*dx2 + dy2*dy2);
        double negCos1 = -dx1 / norm1;
        double negCos2 = -dx2 / norm1;

        if (!float_equal(negCos1,negCos2))
            return negCos1 < negCos2;

        /* If the two angles are the same, return whichever is closer. */
        return norm1 < norm2;
    }
    double orig_x;
    double orig_y;
};

/* An algorithm for finding the convex hull of a set of points in the 2D plane
 * using the Graham Scan.  This algorithm has very good practical runtime
 * (O(n lg n)) and is fairly simple to understand.  The intuition behind the
 * algorithm is to find some point that is known to be on the convex hull of
 * the list of points, then to compute the angles between that point and every
 * other point in the set.  From there, we can sort the points in O(n lg n)
 * time.  The algorithm concludes by marching around these points in sorted
 * order, checking whether each is on the convex hull and adding it if it is.
 *
 * More specifically, the algorithm begins by picking the point with the
 * smallest Y coordinate, which must be on the convex hull.  It then computes
 * the angle made by each other point in the set and the X axis, then sorts
 * these points in ascending order.  Next, we maintain a stack of our current
 * guess of what the convex hull is, which initially is the point with the
 * lowest Y value and the point with the lowest (signed) angle with the x
 * axis.  From there, we iterate across the points in the convex hull
 * expanding our guess.  In particular, let v0 and v1 be the last two points
 * in our convex hull estimate, and let v2 be the next test point.  We then
 * compute the angle between the vector v1 - v0 and the vector v2 - v1.  If
 * this angle is in [pi, 2pi), then we are in a situation like this:
 *
 *                                             /
 *                                v1          /
 *                                *----------*
 *                               /          v0
 *                              /
 *                             /
 *                            *
 *                           v2
 *
 * This means that the point v1 is not actually on the convex hull; it's
 * contained in the hull between v0 and v2.  Of course, depending on what the
 * point in our convex hull estimate is that comes before v0, it's possible
 * that v0 itself isn't on the convex hull either.  In this case, we continue
 * removing nodes from our convex hull estimate until we find that the angle
 * between the last two nodes in the estimate and the new node is in the
 * range [0, pi).  We then update the convex hull by appending this new
 * point.
 *
 * We can now argue the correctness of the algorithm, along with the O(n lg n)
 * runtime.  We can easily see that the algorithm produces a convex polygon,
 * since if we follow the edges from the starting vertex v0 around, each edge
 * turns inward toward the center of the polygon. (A more rigorous proof of
 * this fact exists, but I think this previous one is more intuitive).  To see
 * that the algorithm produces a convex polygon containing all the points in
 * the initial input set, suppose for the sake of contradiction that this
 * isn't true; that there is some node v that isn't in the hull.  It can't be
 * lower than v0, since v0 has the lowest y coordinate, and so it must have
 * been considered by the algorithm at some point.  This means that it was
 * added to the convex hull candidate, and because it wasn't returned it must
 * have been removed at some point.  This means that there was some new point
 * in which the angle between the edge ending at v and the edge containing the
 * new point was in [0, pi).  But by removing this node from consideration,
 * the new edge added between the predecessor of v and the new point has point
 * v in its negative half-space, and so v is in the hull, a contradiction.
 *
 * To argue that the runtime is O(n lg n), we note that we can find the point
 * with the smallest y coordinate in O(n) time, compute the angle each point
 * makes with the x axis in constant time per node (taking O(n) net time), and
 * can then sort them in O(n lg n).  If we can show that the step of growing
 * the hull takes O(n lg n) time, then the overall runtime will be O(n lg n).
 * Initially, it might not seem that this step of the algorithm runs in this
 * time.  Each time we add a node we may have to backtrack all the way through
 * the potentially O(n) nodes on the convex hull, and since we're considering
 * O(n) nodes, this takes O(n^2) time.  However, this analysis is not tight.
 * Notice that once we remove a node from the stack, no future backtracking
 * can ever visit this node again.  This means that in all of the backtracking
 * operations, we can remove at most O(n) nodes from the stack.  Since each
 * backtracking operation takes time linear in the number of nodes removed,
 * this means that the total runtime for all of the backtracking operation is
 * O(n), giving us an overall runtime of O(n lg n).
 */

 template <typename ForwardIterator, typename OutputIterator>
 OutputIterator GrahamScan(ForwardIterator begin, ForwardIterator end,
         OutputIterator out) {

     /* Edge cases - if the range has fewer than three elements, the convex hull
      * is just those points.
      */
     if (size_t(std::distance(begin, end)) < 3)
         return std::copy(begin, end, out);

     /* Locate the element with the smallest y value, breaking ties by choosing
      * coordinates as far to the right (-x) as possible.
      */
     ForwardIterator minY = std::min_element(begin, end, compare_dpoint);

     /* Get an iterator one step past minY; it's the start of the sequence of
      * values that come after it in the input.
      */
     ForwardIterator next = minY; ++next;

     /* We now need to sort the points by their angle with the X axis.  Because
      * we aren't allowed to rearrange the input sequence, we'll make a local
      * copy of the sequence, then will sort that.  We'll leave the lowest point
      * out of the copy so that we don't end up including it in the result.
      */
     DPointVec points;
     points.insert(points.end(), begin, minY); // First portion of the points.
     points.insert(points.end(), next, end);   // Remainder of the points.

     /* Sort by angle with the X axis.  To avoid issues where two adjacent points
      * in the sequence have an 180 degree angle between them, break ties by
      * choosing the point closest to the bottommost point.
      */
     std::sort(points.begin(), points.end(), DPointAngleComparator(minY->lng, minY->lat));

     /* For simplicity, add the minimum point onto the end of the ordering.  This
      * allows us to confirm that the last point we add in the sweep is correct
      * without having to special-case it.
      */
     points.push_back(*minY);

     /* Now, start building up the list of the points in the convex hull.
      * Initially this is the lowest point and the point with the lowest angle,
      * which happens to be the first element of the sorted sequence.
      */
     DPointVec result;
     result.push_back(*minY);
     result.push_back(points[0]);

     /* Now, continuously refine the convex hull until we end up coming back
      * around to the beginning.
      */
     for (size_t i = 1; i < points.size(); ++i) {
         /* Expand the convex hull by factoring in this next point.  This may
          * entail removing some of our previous points, but it always ends by
          * adding this new point.
          */
         while (true) {
             /* Compute two vectors - one consisting of the last two points of the
              * candidate hull, and one consisting of of the last point and the next
              * point in the list.
              */
             double last_x = result[result.size() - 1].lng - result[result.size() - 2].lng;
             double last_y = result[result.size() - 1].lat - result[result.size() - 2].lat;
             double curr_x = points[i].lng - result[result.size() - 1].lng;
             double curr_y = points[i].lat - result[result.size() - 1].lat;

             /* Check whether the angle between these vectors is in the range [0, pi)
              * or [pi, 2*pi).  If it's in the first group, we can add it.  Otherwise
              * we need to remove the last point from the hull and try again.
              *
              * Rather than directly computing the angle between the two vectors, we
              * can instead compute the sine of the angle.  If it's between [0, pi)
              * this will be nonnegative, and if it's between [pi, 2*pi) this would
              * be negative.
              *
              * We can compute the sine of the angle between the vectors by using the
              * 2D cross-product:
              *
              *             |   1   1   1 |
              *   |A x B| = | A.x A.y   0 | = A.x B.y - A.y B.x = |A| |B| sin(theta)
              *             | B.x B.y   0 |
              *
              * Since |A| |B| >= 0, this quantity is positive iff sin(theta) is
              * positive.
              */
             if (last_x * curr_y - last_y * curr_x >= 0) break;

             /* If we're here, it means that this angle was negative and so our last
              * point isn't going to work.  Undo it.
              */
             result.pop_back();
         }

         /* Finally, add the point. */
         result.push_back(points[i]);
     }

     /* At the very end, we now have our convex hull, with the lowest point added
      * twice.  We'll get rid of this point, then return the hull we found.
      */
     result.pop_back();

     /* Move the hull into the output range, then return the endpoint. */
     return std::copy(result.begin(), result.end(), out);
 }

// An algorithm for finding the *median* of a set of points in the 2D plane.
Vector<2> quick_find_median(const std::vector<Vector<2>> &points) {
    std::vector<S2CellId> cells;
    for(auto &pt: points) {
        S2CellId sid = S2CellId::FromLatLng(S2LatLng::FromDegrees(pt[0], pt[1]));
        cells.push_back(sid);
    }
    S2CellId median = cells[cells.size()/2];
    S2LatLng ll = median.ToLatLng();
    Vector<2> result;
    result[0] = ll.lat().degrees();
    result[1] = ll.lng().degrees();
    return result;
}

