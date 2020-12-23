import ConstrainedSystems: recursivecopy!



@testset "Recursive copy" begin

  struct MyStruct{T}
    a :: T
  end

  n = 1000
  a = rand(n,n)
  z = similar(a)

  s1 = MyStruct(a)
  s2 = MyStruct(z)

  recursivecopy!(s2,s1)
  @test s2.a == s1.a
  @test !(s2.a === s1.a)
  @test @allocated(recursivecopy!(s2,s1)) == 16

  ss1 = MyStruct(s1)
  ss2 = MyStruct(MyStruct(z))
  recursivecopy!(ss2,ss1)
  @test ss2.a.a == ss1.a.a
  @test !(ss2.a.a === ss1.a.a)
  @test @allocated(recursivecopy!(s2,s1)) == 16

  t1 = [s1,s1]
  t2 = [MyStruct(z),MyStruct(z)]
  recursivecopy!(t2,t1)
  @test t2[1].a == t1[1].a
  @test t2[2].a == t1[2].a
  @test !(t2[1].a === t1[1].a)
  @test @allocated(recursivecopy!(s2,s1)) == 16


end