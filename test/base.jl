@testset "dot" begin
    
    @test P.dot(zeros(2), ones(2)) ≈ 0.0;
    @test P.dot(ones(3), ones(3)) ≈ 3.0;
    @test P.dot(ones(4), 2*ones(4)) ≈ 8.0;

    S = BitVector([true, false, true]);
    @test P.dot(ones(3), zeros(3), S) ≈ 0.0;
    @test P.dot(ones(3), 2*ones(3), S) ≈ 4.0;
    @test P.dot(ones(3), ones(3), falses(3)) ≈ 0.0;

end

@testset "norm" begin
    
    @test P.norm(zeros(2)) ≈ 0.0;
    @test P.norm(ones(3)) ≈ sqrt(3);
    @test P.norm([3.0, 4.0]) ≈ 5.0;

    S = BitVector([true, false, true]);
    @test P.norm(zeros(3), S) ≈ 0.0;
    @test P.norm(ones(3), S) ≈ sqrt(2);
    @test P.norm([3.0, 3.5, 4.0], S) ≈ 5.0;

end

@testset "project!" begin
    
    x = zeros(4);
    ℓ = [-1.0, 0.5, 0.0, -Inf];
    u = [-0.5, 1.0, 0.0, 0.0];
    P.project!(x, ℓ, u);
    @test x ≈ [-0.5, 0.5, 0.0, 0.0];

end