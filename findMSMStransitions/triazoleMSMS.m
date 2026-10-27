classdef triazoleMSMS
   properties
      name char
      theoreticalMass double
      noOccurences int32 = 0
      collisionEnergy int32 = 0
      retentionTimes (:, 1) double
      foundMass(:, 1) double
      fragmentMass(:, 12) double
      fragmentIntensity(:, 12) double
   end   
end