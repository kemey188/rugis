
CCFLAGS = -g -O2 -fPIC -std=c++11 \
					-DNDEBUG \
					-DARCH_K8 \
					-DS2_USE_EXACTFLOAT \
					-DHASH_NAMESPACE=__gnu_cxx \
					-I. -Igoogle

SOURCES = s2/s1angle.cc \
					s2/s1interval.cc \
					s2/s2.cc \
					s2/s2cap.cc \
					s2/s2cell.cc \
					s2/s2cellid.cc \
					s2/s2cellunion.cc \
					s2/s2edgeindex.cc \
					s2/s2edgeutil.cc \
					s2/s2latlng.cc \
					s2/s2latlngrect.cc \
					s2/s2loop.cc \
					s2/s2pointregion.cc \
					s2/s2polygon.cc \
					s2/s2polygonbuilder.cc \
					s2/s2polyline.cc \
					s2/s2r2rect.cc \
					s2/s2region.cc \
					s2/s2regioncoverer.cc \
					s2/s2regionintersection.cc \
					s2/s2regionunion.cc \
					strings/stringprintf.cc \
					strings/ascii_ctype.cc \
					strings/split.cc \
					strings/strutil.cc \
					base/stringprintf.cc \
					base/logging.cc \
					base/int128.cc \
					base/strtoint.cc \
					util/math/mathlimits.cc \
					util/math/mathutil.cc \
					util/math/exactfloat/exactfloat.cc \
					util/coding/varint.cc \
					util/coding/coder.cc \
					util/hash/hash.cc

all: $(SOURCES)
	g++ -pthread -shared -o libs2.so $(CCFLAGS) $^ -lcrypto
