import mesh
import pic
import Species.species
import numpy

trackers1 = numpy.asarray([0, 4, 15, 26, 136, 900, 901, 1000, 4000, 4001])
ind1 = numpy.asarray([0, 1, 2, 5, 61, 124, 125, 126, 136, 1000])

ind_c = 0
tracker_c = 0
n = 0
while ind_c != len(ind1) and tracker_c != len(trackers1):
    if ind1[ind_c] < trackers1[tracker_c]:
        n += 1
        ind_c += 1
        continue
    elif ind1[ind_c] > trackers1[tracker_c]:
        trackers1[tracker_c] -= n
        tracker_c += 1
        continue
    elif ind1[ind_c] == trackers1[tracker_c]:
        trackers1[tracker_c] = len(ind1)
        n += 1
        tracker_c += 1
        ind_c += 1

if ind_c == len(ind1) and tracker_c < len(trackers1):
    trackers1[tracker_c:] -= n

print(trackers1, tracker_c, ind_c, n)
