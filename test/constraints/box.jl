@testset "project!(x, ::Box)" begin
    
    ℓ = [-1.0, 0.5, 0.0, -Inf];
    u = [-0.5, 1.0, 0.0, 0.0];
    box = Box(ℓ, u);
    
    x = zeros(4);
    P.project!(x, box);
    @test x ≈ [-0.5, 0.5, 0.0, 0.0];
end