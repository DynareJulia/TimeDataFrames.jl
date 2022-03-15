tdf1 = TimeDataFrame(DataFrame(randn(4, 3), :auto), SemesterDate(2000,2))
tdf2 = TimeDataFrame(DataFrame(randn(4, 3), :auto), SemesterDate(2001,1))
tdf3 = TimeDataFrame(DataFrame(randn(4, 3), :auto), SemesterDate(1999,2))
tdf4 = TimeDataFrame(DataFrame(randn(5, 3), :auto), SemesterDate(2002,2))
tdf5 = TimeDataFrame(DataFrame(randn(5, 3), :auto), SemesterDate(2004,2))
tdf6 = TimeDataFrame(DataFrame(randn(5, 3), :auto), SemesterDate(1996,2))

tdf1_orig = copy(tdf1)
tdf2_orig = copy(tdf2)
tdf3_orig = copy(tdf3)
tdf4_orig = copy(tdf4)
tdf5_orig = copy(tdf5)
tdf6_orig = copy(tdf6)
TimeDataFrames.align!(tdf1, tdf1)
@test tdf1 == tdf1_orig

TimeDataFrames.align!(tdf1, tdf2)
@test size(tdf1) == (5, 3)
@test size(tdf2) == (5, 3)
@test periods(tdf1)[1] == SemesterDate(2000,2)
@test periods(tdf1)[end] == SemesterDate(2002,2)
@test periods(tdf2)[1] == SemesterDate(2000,2)
@test periods(tdf2)[end] == SemesterDate(2002,2)
@test tdf1[1, 1] == tdf1_orig[1, 1]
@test tdf1[4, 1] == tdf1_orig[4, 1]
@test ismissing.(tdf1[5, 1])
@test tdf2[2, 1] == tdf2_orig[1, 1]
@test tdf2[5, 1] == tdf2_orig[4, 1]
@test ismissing.(tdf2[1, 1])
