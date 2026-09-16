namespace subdiv{
    namespace curve{
        template<> 
        const std::array<int,2> 
        subdiv::curve::ShapeProp<corlib::LINE>::next = { {1, 0} };

        template<> 
        const std::array<int,2> 
        subdiv::curve::ShapeProp<corlib::LINE>::previous = { {1, 0} };
    }
}
