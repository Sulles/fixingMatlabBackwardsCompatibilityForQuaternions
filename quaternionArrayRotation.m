function [rotated_vector_array] = quaternionArrayRotation(rotation_quaternion_array, vector_array)
	%% Rotate an array of vectors, "vector_array", by an array of rotation quaternions, "rotation_quaternion_array"
	% NOTE: This was tested to be compatible with MATLAB 2019a
	%
	% ARGUMENTS:
	% 	- rotation_quaternion_array: an array of rotation quaternions
	%		+ size is expected to be a 4xN
	%		+ quaternion format is expected to be: q0/r + q1/i + q2/j + q3/k where
	%			r is the "real" component of the quaternion,
	%			i,j,k are the "imaginary" components of the quaternion
	%		+ example of constructing a rotation quaternion array:
	%			num_of_rotations = 10;
	%			rotation_axis_array ("z" axis) = [
	%				zeros(1,num_of_rotations);
	%				zeros(1,num_of_rotations);
	%				linspace(0,1,num_of_rotations)
	%			];
	%			rotation_amount (here 30 degrees) = deg2rad(30);
	%			rotation_quaternion_array = [
	%				cos(rotation_amount./2);
	%				sin(rotation_amount./2) .* rotation_axis_array(1,:);
	%				sin(rotation_amount./2) .* rotation_axis_array(2,:);
	%				sin(rotation_amount./2) .* rotation_axis_array(3,:);
	%			];
	% 	- vector_array: an array of vectors that need to be rotated
	%		+ size is expected to be a 3xN
	%
	% RETURNS:
	% 	- rotated_vector_array: the rotated array of vectors
	%		+ size is expected to be a 3xN
	
	q0 = rotation_quaternion_array(1,:);
	q1 = rotation_quaternion_array(2,:);
	q2 = rotation_quaternion_array(3,:);
	q3 = rotation_quaternion_array(4,:);
	
	% https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
	R1 = [2*(q0.^2 + q1.^2) - 1; 2*(q1.*q2 + q0.*q3); 2*(q1.*q3 - q0.*q2)];
	R2 = [2*(q1.*q2 - q0.*q3); 2*(q0.^2 + q2.^2) - 1; 2*(q2.*q3 + q0.*q1)];
	R3 = [2*(q1.*q3 + q0.*q2); 2*(q2.*q3 - q0.*q1); 2*(q0.^2 + q3.^2) - 1];
	
	rotated_x_component = R1.*vector_array(1,:);
	rotated_y_component = R2.*vector_array(2,:);
	rotated_z_component = R3.*vector_array(3,:);
	
	rotated_vector_array = rotated_x_component + rotated_y_component + rotated_z_component;
end
