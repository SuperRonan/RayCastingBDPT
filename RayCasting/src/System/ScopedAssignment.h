#pragma once

////////////////////////////////////////////////////////
//
// Source: PBRT V3
//
//
////////////////////////////////////////////////////////
template <typename Type>
class ScopedAssignment {

protected:
	Type* target, backup;

public:
	
	ScopedAssignment(Type* target = nullptr, Type value = Type()): 
		target(target) 
	{
		if (target) 
		{
			backup = *target;
			*target = value;
		}
	}
	
	~ScopedAssignment() 
	{
		if (target)* target = backup;
	}

	ScopedAssignment(const ScopedAssignment&) = delete;

	ScopedAssignment& operator=(const ScopedAssignment&) = delete;
	
	ScopedAssignment& operator=(ScopedAssignment&& other) 
	{
		if (target)* target = backup;
		target = other.target;
		backup = other.backup;
		other.target = nullptr;
		return *this;
	}

};